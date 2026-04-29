%% 1. Case Study: Structured bond
clc; close all; clear all;
% Set strikes
strikes = [1.50, 1.75, 2.00, 2.25, 2.50, 3.00, 3.50, 4.00, 5.00, 6.00, 7.00, 8.00, 10.00];

% Set maturities (years)
maturities = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20];
n_strikes   = length(strikes);
n_maturities = length(maturities);


% Init volatility matrix
flat_vols = [
    14.0, 13.0, 12.9, 12.1, 13.3, 13.8, 14.4, 15.0, 17.2, 19.1, 20.2, 21.6, 23.9;
    22.4, 19.7, 17.5, 18.0, 19.2, 20.4, 21.0, 21.4, 22.3, 23.6, 24.9, 26.1, 28.1;
    23.8, 21.7, 20.0, 19.8, 20.3, 20.5, 20.8, 21.4, 22.9, 24.3, 25.6, 26.7, 28.2;
    24.2, 22.4, 20.9, 20.4, 20.4, 20.2, 20.2, 20.5, 21.7, 22.9, 24.0, 25.0, 26.6;
    24.3, 22.6, 21.2, 20.6, 20.4, 19.8, 19.5, 19.6, 20.5, 21.5, 22.6, 23.5, 25.0;
    24.3, 22.7, 21.4, 20.7, 20.2, 19.4, 18.9, 18.8, 19.3, 20.2, 21.2, 22.0, 23.5;
    24.1, 22.6, 21.4, 20.7, 20.1, 19.1, 18.4, 18.1, 18.4, 19.1, 20.0, 20.8, 22.2;
    23.9, 22.5, 21.4, 20.6, 20.0, 18.8, 18.0, 17.6, 17.6, 18.2, 19.0, 19.8, 21.1;
    23.7, 22.4, 21.3, 20.5, 19.8, 18.5, 17.6, 17.1, 17.0, 17.6, 18.3, 19.0, 20.3;
    23.5, 22.2, 21.2, 20.4, 19.6, 18.3, 17.3, 16.8, 16.5, 16.9, 17.6, 18.3, 19.5;
    23.0, 21.7, 20.8, 20.0, 19.3, 17.9, 16.9, 16.2, 15.8, 16.0, 16.5, 17.1, 18.1;
    22.3, 21.2, 20.3, 19.5, 18.7, 17.3, 16.3, 15.5, 15.0, 15.1, 15.5, 16.0, 16.9;
    21.6, 20.4, 19.5, 18.8, 18.0, 16.6, 15.5, 14.7, 14.1, 14.1, 14.5, 15.0, 15.9
];

flat_vols = flat_vols / 100;   % -> convert vols to decimals
% Load and bootstrap the zero rate curve
formatDate = 'dd/mm/yyyy';
[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet);
t0 = dates(1);  % Settlement date: 19/02/2008

% Exercise 1a: LMM Spot Vol Calibration
%
% Build the quarterly caplet schedule, retrieve discount factors and
% forward rates on it, then bootstrap the LMM spot vols on the same
% (strike, maturity) grid as the input flat vols. The first caplet
% (t0 -> 3m) is treated as already fixed at t0, so it contributes only
% an intrinsic-value term to every cap.

%%  SECTION 1: Build the quarterly caplet schedule

% Raw quarterly dates from t0 + 3m up to t0 + 20Y. addtodate is scalar,
% so arrayfun is used to vectorize over the quarter index.
n_quarters = 4 * 20;
raw_dates  = arrayfun(@(q) addtodate(t0, 3*q, 'month'), (1:n_quarters)');

% Adjust each raw date with the Modified Following BD convention.
cap_dates   = adj_modfollow(raw_dates); % first reset 19-May-2008, last 22-Feb-2028
n_cap_dates = length(cap_dates);

% The discount curve is sampled at irregular dates (depos / futures /
% swap pillars), so we resample it onto cap_dates in the next section.

%% SECTION 2: Compute discount factors and forward rates on the caplet schedule
% Discount factors at every caplet reset date.
B_cap = linearRateInterp(dates, discounts, t0, cap_dates); % B(t0, T_i) for each caplet

n_caplets = length(cap_dates); % number of caplet payment dates

% Pre-allocation (overwritten below by the vectorized assignments).
fwd_rates  = zeros(n_caplets, 1); % forward rate of each caplet
delta_fwd  = zeros(n_caplets, 1); % Act/360 year fraction of each caplet period
tau_expiry = zeros(n_caplets, 1); % Act/365 time from t0 to caplet reset

% Reset / payment dates (column vectors).
T_reset   = cap_dates(1:end-1);
T_payment = cap_dates(2:end);

% Year fractions via MATLAB built-ins (basis 2 = Act/360, basis 3 = Act/365).
delta_fwd  = yearfrac(T_reset, T_payment, 2);
tau_expiry = yearfrac(t0,      T_reset,   3);   % time from t0 to each caplet reset

% Discount factors at the reset and payment date of each caplet.
B_Ti  = B_cap(1:end-1);
B_Ti1 = B_cap(2:end);

% Quarterly forward Euribor rates implied by the discount factors.
fwd_rates = (B_Ti ./ B_Ti1 - 1) ./ delta_fwd;

%% SECTION 3: Map cap maturities to caplet indices

% Target dates = t0 + n years, then Modified Following adjustment.
raw_targets  = arrayfun(@(y) addtodate(t0, y, 'year'), maturities(:));
target_dates = adj_modfollow(raw_targets);

% Position of each target date inside cap_dates (vectorized search).
% After the -1 shift below, cap_maturity_idx(j) is the index of the
% LAST caplet that belongs to cap j.
[found, cap_maturity_idx] = ismember(target_dates, cap_dates);
if ~all(found)
    error('Cap maturities not found in cap_dates: %s', ...
          mat2str(maturities(~found)));
end
cap_maturity_idx = cap_maturity_idx -1;% index of the last caplet of each cap
%% SECTION 4: Bootstrap LMM spot vols
%  Per strike, iterate over the maturity buckets:
%    1. Re-price the cap at T_j with the market flat vol.
%    2. Subtract the price of caplets already pinned in earlier buckets.
%    3. Solve via fzero for the right-edge sigma_beta so the new caplets
%       (with linearly-interpolated vols between sigma_alpha and
%       sigma_beta) reproduce that residual.
%    4. Store the interpolated spot vols of the bucket.
%
%  The 1Y bucket is pre-filled with the 1Y flat vol before entering the
%  loop, so the bucket loop starts at j = 2.
%
%  OUTPUT: spot_vols(i, k) = spot vol of caplet i for strike k.

% Pre-allocate the spot vol matrix.
% Rows = caplet index, Cols = strike index.
n_total_caplets = cap_maturity_idx(end);
spot_vols = zeros(n_total_caplets, n_strikes);

for j = 1:cap_maturity_idx(1)
spot_vols(j, :) = flat_vols(1,:);
end

% Loop over strikes
for k = 1:n_strikes
    K = strikes(k) / 100; % strike in decimal

    % Initial state for this strike: every caplet in the 1Y bucket shares
    % the 1Y flat vol, so we read it from the pre-filled spot_vols.
    sigma_alpha = spot_vols(1, k);

    idx_all = (1:cap_maturity_idx(1))';
    % B_cap(idx_all+1) gives the discount factor at the payment date of
    % each caplet in the 1Y bucket.
    cap_price_already_priced = sum( caplet_black_LMM( ...
                                    fwd_rates(idx_all), K, delta_fwd(idx_all), ...
                                    B_cap(idx_all+1),  tau_expiry(idx_all), sigma_alpha) );

    for j = 2:n_maturities
        Sigma_beta = flat_vols(j, k);           % market flat vol for cap j at strike k
        i_alpha = cap_maturity_idx(j-1)+1;      % first new caplet in this bucket
        i_beta  = cap_maturity_idx(j);          % last caplet of cap j

        % --- Step 1: market cap price at T_j (vectorized over caplets) ---
        idx_all = (1:i_beta)'; % every caplet that composes cap j
        cap_price_flat = sum( caplet_black_LMM( ...
            fwd_rates(idx_all), K, delta_fwd(idx_all), ...
            B_cap(idx_all+1),  tau_expiry(idx_all), Sigma_beta) );

        % --- Step 2: price of NEW caplets only ---
        delta_cap_price = cap_price_flat - cap_price_already_priced;

        % --- Step 3: bucket boundaries and root finding ---
        % T_alpha = reset of the last caplet of the previous cap,
        % T_beta  = reset of the last caplet of the current cap.
        T_alpha = T_reset(i_alpha-1);
        T_beta = T_reset(i_beta);

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

        % --- Step 4: store spot vols for this bucket (vectorized) ---
        idx_b = (i_alpha:i_beta)';
        w = (T_reset(idx_b)- T_alpha) ./ (T_beta - T_alpha);
        spot_vols(idx_b, k) = sigma_alpha + w .* (sigma_beta - sigma_alpha);
    

        % --- Step 5: update running totals for the next bucket ---
        % After fzero converges, the cap price already covered up to
        % i_beta equals cap_price_flat, so no recomputation is needed.
        sigma_alpha              = sigma_beta;
        cap_price_already_priced = cap_price_flat;
    end
end
%% SECTION 5: Display results
%  Print the spot vols at cap maturity dates so the output sits on the
%  same grid as the input flat vols for direct comparison.

fprintf('\n=== LMM Spot Vols (at cap maturity dates) [%%] ===\n');
fprintf('Maturity |');
for k = 1:n_strikes
    fprintf(' K=%.2f%%', strikes(k));
end
fprintf('\n');

for j = 1:n_maturities
    fprintf('  %3dy   |', maturities(j));
    i_beta = cap_maturity_idx(j);
    for k = 1:n_strikes
        fprintf('  %6.2f ', spot_vols(i_beta, k) * 100);
    end
    fprintf('\n');
end