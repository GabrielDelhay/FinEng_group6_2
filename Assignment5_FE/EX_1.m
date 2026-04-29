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
%  GOAL: Extract LMM spot vols from market flat cap vols
%
%  THEORY:
%    - A Cap(T,K) = sum of Caplets, each priced with Black:
%        caplet_i = B(0,T_{i+1}) * delta_i * [L_i*N(d1) - K*N(d2)]
%      where d1,d2 use the SPOT vol sigma_i (not flat vol)
%    - Flat vol Sigma_n: same vol for ALL caplets in cap -> 
%        Cap(T_n, K) = sum_i caplet(i, Sigma_n)   [market quote]
%    - Spot vol sigma_i: TRUE vol of forward rate i ->
%        Cap(T_n, K) = sum_i caplet(i, sigma_i)   [model decomp]
%    - Bootstrap: for each new maturity T_n, the price of the
%      NEW caplets = Cap(T_n) - Cap(T_{n-1}), so we solve for
%      the spot vols in the new bucket via the linear constraint
%      (slide 23):
%        sigma_i = sigma(T_alpha) + (T_i - T_alpha)/(T_beta - T_alpha)
%                  * (sigma(T_beta) - sigma(T_alpha))
%      i.e. spot vols are LINEARLY interpolated within each bucket
%
% the caplet from T_0 to T_1 (first 3m) is assumed already fixed at t_0, so caps start from caplet 2.

%%  SECTION 1: Build the quarterly caplet schedule

% Quarterly raw dates (3m, 6m, ..., 20Y) — addtodate is scalar, so arrayfun
n_quarters = 4 * 20; % don't consider the last couplet
raw_dates  = arrayfun(@(q) addtodate(t0, 3*q, 'month'), (1:n_quarters)');

% Apply Modified Following BD adjustment in one shot (vectorized)
cap_dates   = adj_modfollow(raw_dates); % from 19-May-2008 to 22-Feb-2028: so is the reset date for each caplet
n_cap_dates = length(cap_dates);

% The bootstrapped dates are NOT purely quarterly. We therefore
% interpolate discount factors from the bootstrap curve onto
% the quarterly caplet schedule defined above.
% log-linear interpolation is used (standard for discount factors).

%% SECTION 2: Compute discount factors and forward rates on the caplet schedule
% Interpolate discount factors at caplet dates
B_cap = linearRateInterp(dates, discounts, t0, cap_dates); %discounts in the caplets reset dates

n_caplets = length(cap_dates); % number of potential payment dates

% Pre-allocate forward rates and day count fractions
fwd_rates = zeros(n_caplets, 1); % L_i = forward rate for caplet i
delta_fwd = zeros(n_caplets, 1); % Act/360 year fraction for coupon
tau_expiry = zeros(n_caplets, 1); % Act/365 year fraction for Black

% I DAYCOUNT DEVO RIGUARDARLI

% Reset / payment dates (column vectors)
T_reset   = cap_dates(1:end-1);
T_payment = cap_dates(2:end);

% Year fractions via MATLAB built-ins (basis 2 = Act/360, basis 3 = Act/365)
delta_fwd  = yearfrac(T_reset, T_payment, 2);  
tau_expiry = yearfrac(t0,      T_reset,   3);   % yearfrac between t0 and Libor reset time

% Discount factors a T_i e T_{i+1}
B_Ti  = B_cap(1:end-1);
B_Ti1 = B_cap(2:end);

% Forward Euribor 3m: L_i = (1/delta_i) * (B(0,T_i)/B(0,T_{i+1}) - 1)
fwd_rates = (B_Ti ./ B_Ti1 - 1) ./ delta_fwd;

%% SECTION 3: Map cap maturities to caplet indices
%
%  For each cap maturity (1Y, 2Y, ..., 20Y) find the caplet index that
%  pays at that date, applying the Modified Following BD Convention.
%  (First caplet excluded: cap starts from i_start = 2)

% Target dates = t0 + n years, then Modified Following adjustment
raw_targets  = arrayfun(@(y) addtodate(t0, y, 'year'), maturities(:));
target_dates = adj_modfollow(raw_targets);

% Corresponding indices inside cap_dates (vectorized)
% cap_maturity_indx signals the first caplet of the new cap which has to
% be computed
[found, cap_maturity_idx] = ismember(target_dates, cap_dates);
if ~all(found)
    error('Cap maturities not found in cap_dates: %s', ...
          mat2str(maturities(~found)));
end
cap_maturity_idx = cap_maturity_idx -1;% representing the index of the last couplet in each cap
%% SECTION 4: Bootstrap LMM spot vols
%  ALGORITHM:
%  For each strike k:
%    For each maturity bucket [T_{alpha}, T_{beta}]:
%      1. Compute market cap price Cap_mkt(T_beta, K) using flat vol
%      2. Subtract already-priced caplets from previous buckets
%         -> get Delta_C = price of NEW caplets in this bucket
%      3. Within the bucket, spot vols are LINEARLY interpolated:
%           sigma_i = sigma_alpha + (T_i - T_alpha)/(T_beta - T_alpha)
%                     * (sigma_beta - sigma_alpha)
%         where sigma_alpha is known (last bucket's endpoint)
%         and sigma_beta is the unknown to solve for
%      4. Solve for sigma_beta such that 
%           sum_{i in bucket} caplet(i, sigma_i) = Delta_C
%         using fzero (1D root finding)
%
%  OUTPUT: spot_vols(i, k) = spot vol of caplet i for strike k
%          (same grid as flat vols: 13 strikes x n_caplets)

% Pre-allocate spot vol matrix
% Rows = caplet index (from i_start to last cap maturity)
% Cols = strike index
n_total_caplets = cap_maturity_idx(end);
spot_vols = zeros(n_total_caplets, n_strikes);

for j = 1:cap_maturity_idx(1)
spot_vols(j, :) = flat_vols(1,:);
end

% Loop over strikes
for k = 1:n_strikes
    K = strikes(k) / 100; % strike in decimal

    % Initial state for this strike
    sigma_alpha = spot_vols(1, k);             

    idx_all = (1:cap_maturity_idx(1))';
    % we use  B_cap(idx_all+1) in order to obtain DF in payment date of each caplet
    cap_price_already_priced = sum( caplet_black_LMM( ...
                                    fwd_rates(idx_all), K, delta_fwd(idx_all), ...
                                    B_cap(idx_all+1),  tau_expiry(idx_all), sigma_alpha) );

    for j = 2:n_maturities
        Sigma_beta = flat_vols(j, k);           % flat vol for this (mat, strike)
        i_alpha = cap_maturity_idx(j-1)+1;
        i_beta  = cap_maturity_idx(j);        % last caplet for cap j

        % --- Step 1: market cap price at T_j (vectorized over caplets) ---
        idx_all = (1:i_beta)'; % caplets composing the new cap
        cap_price_flat = sum( caplet_black_LMM( ...
            fwd_rates(idx_all), K, delta_fwd(idx_all), ...
            B_cap(idx_all+1),  tau_expiry(idx_all), Sigma_beta) );

        % --- Step 2: price of NEW caplets only ---
        delta_cap_price = cap_price_flat - cap_price_already_priced;

        % --- Step 3: bucket boundaries and root finding ---
     

%occhio a alscaire vuoto
%datetime(T_payment, 'ConvertFrom', 'datenum')
        %those are the initial and final times of our CAP
        T_alpha = T_reset(i_alpha-1); % taking the index of the last caplet in the previous cap
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
    

        % --- Step 5: update running totals (NO recomputation needed) ---
        % After fzero converges, sum of caplets in (i_start..i_beta_new) with
        % the stored spot vols is exactly cap_price_flat by construction.
        sigma_alpha              = sigma_beta;
        cap_price_already_priced = cap_price_flat;
    end
end
%% SECTION 5: Display results
%  Show spot vols at cap maturity dates (same grid as input) for easy comparison with flat vols

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