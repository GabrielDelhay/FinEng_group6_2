%% 1. Case Study: Structured bond
clc; close all; clear all;
% Set strikes
strikes = [1.50, 1.75, 2.00, 2.25, 2.50, 3.00, 3.50, 4.00, 5.00, 6.00, 7.00, 8.00, 10.00];

% Set maturities (years)
maturities = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20];
n_strikes   = length(strikes);
n_maturities = length(maturities);
bp = 0.0001;

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
formatDate = 'dd/mm/yyyy';
[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet);
t0 = dates(1);  % Settlement date: 19/02/2008


%% EXERCISE 1.a
% Load and bootstrap the zero rate curve
[spot_vols, cap_dates, B_cap, fwd_rates, delta_fwd, tau_expiry, T_reset, cap_maturity_idx] = ...
    lmm_spot_vols(flat_vols, strikes, maturities, dates, discounts, t0);
%% Display results
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

%% EXERCISE 1.b and 1.g
%% Bond parameters from the termsheet
N        = 50e6;       % notional, 50 MIO EUR
spread   = 0.0200;       % Party A pays Euribor 3m + 2.00%

bond.first_cpn = 0.04;
bond.spread1   = 0.0100;   bond.K1 = 0.0420;   bond.c_dig1 = 0.0070;  % Regime 1
bond.spread2   = 0.0120;   bond.K2 = 0.0470;   bond.c_dig2 = 0.0100;  % Regime 2 
bond.spread3   = 0.0130;   bond.K3 = 0.0540;   bond.c_dig3 = 0.0110;  % Regime 3

%% Calculate pricing (Flat Black and Smile-corrected)
[NPV_A, NPV_B_flat, X_flat, NPV_B_smile, X_smile] = ...
    price_structured_bond(N, spread, bond, B_cap, delta_fwd, tau_expiry, ...
                          fwd_rates, cap_maturity_idx, spot_vols, strikes);
%% Output
fprintf('\n=== Pricing of the structured bond ===\n');
fprintf('NPV Party A (floating)         : %14.2f EUR\n',  NPV_A);
fprintf('\n--- Point 1.b: flat Black ---\n');
fprintf('NPV Party B (coupons)          : %14.2f EUR\n',  NPV_B_flat);
fprintf('Upfront X*N                    : %14.2f EUR\n',  X_flat*N);
fprintf('Upfront X                      : %10.4f %%\n',   X_flat*100);
fprintf('\n--- Point 1.g: smile-corrected ---\n');
fprintf('NPV Party B (coupons)          : %14.2f EUR\n',  NPV_B_smile);
fprintf('Upfront X*N                    : %14.2f EUR\n',  X_smile*N);
fprintf('Upfront X                      : %10.4f %%\n',   X_smile*100);
fprintf('\n--- Digital risk size ---\n');
fprintf('Delta upfront X                : %10.4f bp\n',   (X_smile - X_flat)*1e4);

%% EXERCISE 1.c
BPV = 0.0001;
n_depos = 3;
n_futures = 7;
n_swaps = 50;
n_total = n_depos + n_futures + (n_swaps - 1);

% Pre-allocate price vector
price = zeros(n_total, 1);

% Compute initial true_price (Baseline) in EUR
[~, ~, X_flat_base, ~, ~] = price_structured_bond(N, spread, bond, B_cap, delta_fwd, tau_expiry, fwd_rates, cap_maturity_idx, spot_vols, strikes);
true_price = X_flat_base * N;

idx = 1; 

% 1. Bump Depos
for i = 1:n_depos
    ratesSet.depos(i, :) = ratesSet.depos(i, :) + BPV;  % Bump
    [dates, disc_bump, ~] = bootstrap(datesSet, ratesSet);
    
    % Recompute curve-dependent variables (forwards, spot vols)
    [spot_bump, ~, B_bump, fwd_bump, delta_bump, tau_bump, ~, idx_bump] = ...
        lmm_spot_vols(flat_vols, strikes, maturities, dates, disc_bump, t0);
    
    % Recompute price and extract new Upfront
    [~, ~, X_flat_bump, ~, ~] = price_structured_bond(N, spread, bond, B_bump, delta_bump, tau_bump, fwd_bump, idx_bump, spot_bump, strikes);
    price(idx) = X_flat_bump * N;
    
    ratesSet.depos(i, :) = ratesSet.depos(i, :) - BPV;  % Restore
    idx = idx + 1;
end

% 2. Bump Futures
for i = 1:n_futures
    ratesSet.futures(i, :) = ratesSet.futures(i, :) + BPV; 
    [dates, disc_bump, ~] = bootstrap(datesSet, ratesSet);
    
    [spot_bump, ~, B_bump, fwd_bump, delta_bump, tau_bump, ~, idx_bump] = ...
        lmm_spot_vols(flat_vols, strikes, maturities, dates, disc_bump, t0);
        
    [~, ~, X_flat_bump, ~, ~] = price_structured_bond(N, spread, bond, B_bump, delta_bump, tau_bump, fwd_bump, idx_bump, spot_bump, strikes);
    price(idx) = X_flat_bump * N;
    
    ratesSet.futures(i, :) = ratesSet.futures(i, :) - BPV;
    idx = idx + 1;
end

% 3. Bump Swaps (skip the first)
for i = 2:n_swaps
    ratesSet.swaps(i, :) = ratesSet.swaps(i, :) + BPV; 
    [dates, disc_bump, ~] = bootstrap(datesSet, ratesSet);
    
    [spot_bump, ~, B_bump, fwd_bump, delta_bump, tau_bump, ~, idx_bump] = ...
        lmm_spot_vols(flat_vols, strikes, maturities, dates, disc_bump, t0);
        
    [~, ~, X_flat_bump, ~, ~] = price_structured_bond(N, spread, bond, B_bump, delta_bump, tau_bump, fwd_bump, idx_bump, spot_bump, strikes);
    price(idx) = X_flat_bump * N;
    
    ratesSet.swaps(i, :) = ratesSet.swaps(i, :) - BPV; 
    idx = idx + 1;
end

% Vectorized delta calculation
delta = price - true_price;

%% Display results
labels = [strcat("Depo_", string(1:n_depos)), strcat("Future_", string(1:n_futures)), strcat("Swap_", string(2:15))].';
T_delta = table(labels, delta(1:length(labels)), 'VariableNames', {'Pillar', 'PV01_EUR'});
fprintf('\n=== Delta Bucket Analysis (PV01) ===\n');
disp(T_delta);
fprintf('Total Delta (Sum of Buckets): %.2f EUR\n', sum(delta));
%% EXERCISE 1.d
% Total vega
dVol = 0.0001;
flat_vols_up = flat_vols + dVol;
[spot_vols_up, ~, ~, ~, ~, ~, ~, ~] = ...
lmm_spot_vols(flat_vols_up, strikes, maturities, dates, discounts, t0);
[~, ~, X_flat_up_vega, ~, ~] = ...
    price_structured_bond(N, spread, bond, B_cap, delta_fwd, tau_expiry, ...
                          fwd_rates, cap_maturity_idx, spot_vols_up, strikes);
flat_vols_down = flat_vols - dVol;
[spot_vols_down, ~, ~, ~, ~, ~, ~, ~] = ...
lmm_spot_vols(flat_vols_down, strikes, maturities, dates, discounts, t0);
[~, ~, X_flat_down_vega, ~, ~] = ...
    price_structured_bond(N, spread, bond, B_cap, delta_fwd, tau_expiry, ...
                          fwd_rates, cap_maturity_idx, spot_vols_down, strikes);
vega_total_b = (X_flat_up_vega - X_flat_down_vega) * N * 100/2; 
fprintf('Total Vega (+1%% of flat vol): %.2f EUR\n', vega_total_b);

%% EXERCISE 1.e
% Coarse-grained bucket deltas (0-2y, 2-6y, 6-10y)
% Bucket 1: 0-2y  -> Depos (1:3) + Futures (1:7) + Swap_2 (swap idx 2)
% Bucket 2: 2-6y  -> Swaps from ~2y to 6y (swap idx 3:6)
% Bucket 3: 6-10y -> Swaps from ~6y to 10y (swap idx 7:10)
[~, ~, X_flat_base, ~, ~] = price_structured_bond(N, spread, bond, B_cap, delta_fwd, tau_expiry, fwd_rates, cap_maturity_idx, spot_vols, strikes);
true_price = X_flat_base * N;
bucket_delta = coarse_bucket_delta(ratesSet, datesSet, dates, flat_vols, strikes, maturities, t0, N, spread, bond, true_price, BPV, n_depos, n_futures);

fprintf(['\n=== Coarse-Grained Bucket Deltas ===\n', 'Bucket 0-2y  : %10.2f EUR\n', 'Bucket 2-6y  : %10.2f EUR\n', ...
         'Bucket 6-10y : %10.2f EUR\n','Total        : %10.2f EUR\n'], bucket_delta(1), bucket_delta(2), bucket_delta(3), sum(bucket_delta));

% MISSING HEDGING WITH SWAPS

%% EXERCISE 1.f
bucket_vega = coarse_bucket_vega(flat_vols, strikes, maturities, dates, discounts, ...
    t0, N, spread, bond, B_cap, delta_fwd, tau_expiry, ...
    fwd_rates, cap_maturity_idx, spot_vols, true_price, dVol);

fprintf('\n=== Coarse-Grained Bucket Vegas ===\n');
fprintf('Bucket 0-6y  : %10.2f EUR\n', bucket_vega(1));
fprintf('Bucket 6-10y : %10.2f EUR\n', bucket_vega(2));
fprintf('Total        : %10.2f EUR\n', sum(bucket_vega));

% MISSING HEDGING WITH CAPS