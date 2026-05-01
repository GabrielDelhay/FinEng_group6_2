function [spot_vols, cap_dates, B_cap, fwd_rates, delta_fwd, ...
          tau_expiry, T_reset, cap_maturity_idx] = ...
          lmm_spot_vols(flat_vols, strikes, maturities, ...
                                  dates, discounts, t0)
% LMM_BOOTSTRAP_SPOT_VOLS  Bootstrap LMM spot vols from market flat cap vols.
%
%   Internally performs:
%     1. Build quarterly caplet schedule (Mod Following)
%     2. Discount factors and forward rates on the schedule
%     3. Map cap maturities to caplet indices
%     4. Bucket-by-bucket spot vol bootstrap (linear interp, fzero)
%
%   INPUT
%     flat_vols   [n_maturities x n_strikes]  market flat vols (decimals)
%     strikes     [1 x n_strikes]             cap strikes (percent, e.g. 2.0)
%     maturities  [1 x n_maturities]          cap maturities in years
%     dates       [n_curve x 1]               curve pillar dates (datenum)
%     discounts   [n_curve x 1]               discount factors at pillars
%     t0          scalar                      settlement date (datenum)
%
%   OUTPUT
%     spot_vols        [n_total_caplets x n_strikes]  LMM spot vols (decimals)
%     cap_dates        [n_cap_dates x 1]              full caplet schedule (datenum)
%     B_cap            [n_cap_dates x 1]              discount factors on schedule
%     fwd_rates        [n_caplets x 1]                quarterly forward Euribor rates
%     delta_fwd        [n_caplets x 1]                Act/360 caplet year fractions
%     tau_expiry       [n_caplets x 1]                Act/365 time from t0 to reset
%     T_reset          [n_caplets x 1]                caplet reset dates (datenum)
%     cap_maturity_idx [n_maturities x 1]             index of last caplet for each cap

n_strikes    = length(strikes);
n_maturities = length(maturities);

%% SECTION 1: Build the quarterly caplet schedule
% Raw quarterly dates from t0+3m up to t0+20Y, then Mod Following adjustment.
n_quarters = 4 * 20;
raw_dates  = arrayfun(@(q) addtodate(t0, 3*q, 'month'), (1:n_quarters)');
cap_dates  = busdate(raw_dates, 'modifiedfollow');

%% SECTION 2: Discount factors and forward rates on the caplet schedule
B_cap = linearRateInterp(dates, discounts, t0, cap_dates);

T_reset   = cap_dates(1:end-1);
T_payment = cap_dates(2:end);

delta_fwd  = yearfrac(T_reset, T_payment, 2);  % Act/360
tau_expiry = yearfrac(t0,      T_reset,   3);  % Act/365 time to each caplet reset

B_Ti  = B_cap(1:end-1);
B_Ti1 = B_cap(2:end);
fwd_rates = (B_Ti ./ B_Ti1 - 1) ./ delta_fwd;

%% SECTION 3: Map cap maturities to caplet indices
raw_targets  = arrayfun(@(y) addtodate(t0, y, 'year'), maturities(:));
target_dates = busdate(raw_targets, 'modifiedfollow');

[found, cap_maturity_idx] = ismember(target_dates, cap_dates);
if ~all(found)
    error('Cap maturities not found in cap_dates: %s', mat2str(maturities(~found)));
end
cap_maturity_idx = cap_maturity_idx - 1;  % index of last caplet of each cap

%% SECTION 4: Bootstrap LMM spot vols
% Per strike, iterate over maturity buckets:
%   1. Re-price cap at T_j with market flat vol.
%   2. Subtract price of caplets already pinned in earlier buckets.
%   3. Solve via fzero for sigma_beta s.t. new caplets (with linearly
%      interpolated vols between sigma_alpha and sigma_beta) match residual.
%   4. Store interpolated spot vols for the bucket.
% The 1Y bucket is pre-filled before the loop (all caplets share the 1Y flat vol).

n_total_caplets = cap_maturity_idx(end);
spot_vols = zeros(n_total_caplets, n_strikes);

% Pre-fill 1Y bucket
for j = 1:cap_maturity_idx(1)
    spot_vols(j, :) = flat_vols(1, :);
end

for k = 1:n_strikes
    K = strikes(k) / 100;

    sigma_alpha = spot_vols(1, k);

    idx_all = (1:cap_maturity_idx(1))';
    cap_price_already_priced = sum(caplet_black_LMM( ...
        fwd_rates(idx_all), K, delta_fwd(idx_all), ...
        B_cap(idx_all+1),   tau_expiry(idx_all), sigma_alpha));

    for j = 2:n_maturities
        Sigma_beta = flat_vols(j, k);
        i_alpha = cap_maturity_idx(j-1) + 1;  % first new caplet in this bucket
        i_beta  = cap_maturity_idx(j);         % last caplet of cap j

        % Step 1: market cap price at T_j
        idx_all = (1:i_beta)';
        cap_price_flat = sum(caplet_black_LMM( ...
            fwd_rates(idx_all), K, delta_fwd(idx_all), ...
            B_cap(idx_all+1),   tau_expiry(idx_all), Sigma_beta));

        % Step 2: residual price for the new caplets only
        delta_cap_price = cap_price_flat - cap_price_already_priced;

        % Step 3: root-find sigma_beta via fzero
        T_alpha = T_reset(i_alpha - 1);
        T_beta  = T_reset(i_beta);

        f = @(sb) sum_caplets_linear_vol(i_alpha, i_beta, ...
                    fwd_rates, K, delta_fwd, B_cap, T_reset, ...
                    sigma_alpha, sb, T_alpha, T_beta, tau_expiry) ...
                  - delta_cap_price;

        try
            sigma_beta = fzero(f, [1e-6, 2.0]);
        catch
            warning('fzero failed for maturity %dy, strike %.2f%%. Using sigma_alpha.', ...
                    maturities(j), strikes(k));
            sigma_beta = sigma_alpha;
        end

        % Step 4: store linearly interpolated spot vols for this bucket
        idx_b = (i_alpha:i_beta)';
        w = (T_reset(idx_b) - T_alpha) ./ (T_beta - T_alpha);
        spot_vols(idx_b, k) = sigma_alpha + w .* (sigma_beta - sigma_alpha);

        % Step 5: update state for next bucket
        sigma_alpha              = sigma_beta;
        cap_price_already_priced = cap_price_flat;
    end
end

end 