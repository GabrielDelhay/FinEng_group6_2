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
                          fwd_rates, cap_maturity_idx, spot_vols, strikes, cap_dates, t0);
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
n_depos   = 3;
n_futures = 7;
n_swaps   = 15;
n_total   = n_depos + n_futures + (n_swaps - 1);


true_price = X_flat * N;
% Calculate all Delta Buckets using the new function
[delta, T_delta] = calc_delta_buckets(N, spread, bond, datesSet, ratesSet, flat_vols, strikes, maturities, t0, true_price,n_depos, n_futures, n_swaps, n_total, BPV, cap_dates);

%% Display results
fprintf('\n=== Delta Bucket Analysis (PV01) ===\n');
disp(T_delta);
fprintf('Total Delta (Sum of Buckets): %.2f EUR\n', sum(delta));
%% EXERCISE 1.d
% Total vega
dVol = 0.01;
flat_vols_up = flat_vols + dVol;
[spot_vols_up, ~, ~, ~, ~, ~, ~, ~] = ...
lmm_spot_vols(flat_vols_up, strikes, maturities, dates, discounts, t0);
[~, ~, X_flat_up_vega, ~, ~] = ...
    price_structured_bond(N, spread, bond, B_cap, delta_fwd, tau_expiry, ...
                          fwd_rates, cap_maturity_idx, spot_vols_up, strikes,cap_dates, t0);
flat_vols_down = flat_vols - dVol;
[spot_vols_down, ~, ~, ~, ~, ~, ~, ~] = ...
lmm_spot_vols(flat_vols_down, strikes, maturities, dates, discounts, t0);
[~, ~, X_flat_down_vega, ~, ~] = ...
    price_structured_bond(N, spread, bond, B_cap, delta_fwd, tau_expiry, ...
                          fwd_rates, cap_maturity_idx, spot_vols_down, strikes, cap_dates, t0);
vega_total_b = (X_flat_up_vega - X_flat_down_vega) * N /(2); 
fprintf('Total Vega (+1%% of flat vol): %.2f EUR\n', vega_total_b);

%% EXERCISE 1.e
% Coarse-grained bucket deltas (0-2y, 2-6y, 6-10y)
% Bucket 1: 0-2y  -> Depos (1:3) + Futures (1:7) + Swap_2 (swap idx 2)
% Bucket 2: 2-6y  -> Swaps from ~2y to 6y (swap idx 3:6)
% Bucket 3: 6-10y -> Swaps from ~6y to 10y (swap idx 7:10)


[bucket_delta, bucket_DF_bump] = coarse_bucket_delta(ratesSet, datesSet, dates, flat_vols, strikes, maturities, t0, N, spread, bond, true_price, BPV, n_depos, n_futures, cap_dates);
fprintf(['\n=== Coarse-Grained Bucket Deltas ===\n', 'Bucket 0-2y  : %10.2f EUR\n', 'Bucket 2-6y  : %10.2f EUR\n', ...
         'Bucket 6-10y : %10.2f EUR\n','Total        : %10.2f EUR\n'], bucket_delta(1), bucket_delta(2), bucket_delta(3), sum(bucket_delta));

T_expiry = [datesSet.swaps(2); datesSet.swaps(6); datesSet.swaps(10)];
n_buckets = 3;

delta_swaps = zeros(n_buckets, n_buckets);
for j = 1:n_buckets                 % loop over hedging swaps
    for i = 1:n_buckets             % loop over coarse buckets
        delta_swaps(i, j) = calc_swap_delta_exact(t0, T_expiry(j), dates, discounts, bucket_DF_bump(:, i));
    end
end

% Matlab will recon it's a triangular back-substitution: longest swap first
N_hedge_delta = delta_swaps \ (-bucket_delta);       
%% EXERCISE 1.f
bucket_vega = coarse_bucket_vega(flat_vols, strikes, maturities, dates, discounts, ...
                                           t0, N, spread, bond, B_cap, delta_fwd, tau_expiry, ...
                                           fwd_rates, cap_maturity_idx, dVol,cap_dates);

fprintf('\n=== Coarse-Grained Bucket Vegas ===\n');
fprintf('Bucket 0-6y  : %.2f EUR\n', bucket_vega(1));
fprintf('Bucket 6-10y : %.2f EUR\n', bucket_vega(2));
fprintf('Total        : %.2f EUR\n', sum(bucket_vega));

% Caps Vega Matrix
vega_caps = calc_hedge_caps_vega(flat_vols, strikes, maturities, dates, discounts, t0, dVol, N);
fprintf('\n=== Caps Vega Matrix ===\n');
disp(array2table(vega_caps, 'VariableNames', {'Cap_6y', 'Cap_10y'}, 'RowNames', {'Bucket_0_6y', 'Bucket_6_10y'}));

% triangular system
N_hedge_vega = vega_caps \ (-bucket_vega);

fprintf('\n=== Vega Hedging Notionals (EUR) ===\n');
fprintf('Notional Cap 6y  : %.2f\n', N_hedge_vega(1));
fprintf('Notional Cap 10y : %.2f\n', N_hedge_vega(2));

%% Case study 2 : Exotic Cap - Calibration & Pricing under Bond Market Model (BMM)
fprintf('\n=== Exercice 2 ===\n');

%% 1. Market data
N          = 16;
resetDates = buildResetDates(t0, N);
delta      = yearfrac(resetDates(1:end-1), resetDates(2:end), 3);
B0_T       = linearRateInterp(dates, discounts, t0, resetDates);
L0         = (B0_T(1:end-1) ./ B0_T(2:end) - 1) ./ delta;

%% 2. Flat ATM vols from market table
maturitiesYears = [1, 2, 3, 4];
strikes_table   = [1.50,1.75,2.00,2.25,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,10.0] / 100;
vols_table = [
    0.140,0.130,0.129,0.121,0.133,0.138,0.144,0.150,0.172,0.191,0.202,0.216,0.239;
    0.224,0.197,0.175,0.180,0.192,0.204,0.210,0.214,0.223,0.236,0.249,0.261,0.281;
    0.238,0.217,0.200,0.198,0.203,0.205,0.208,0.214,0.229,0.243,0.256,0.267,0.282;
    0.242,0.224,0.209,0.204,0.204,0.202,0.202,0.205,0.217,0.229,0.240,0.250,0.266;
];
flatVolsATM_annual = getATMFlatVols(B0_T, delta, maturitiesYears, strikes_table, vols_table);

%% 3. BMM spot vol calibration
nu = spotvolbootstrap(maturitiesYears, flatVolsATM_annual, t0, resetDates, B0_T, delta, L0);

%% 4. Correlation matrix rho_{ij} = exp(-lambda * |T_i - T_j|)
lambda = 0.1;
tGrid  = yearfrac(t0, resetDates(2:end), 3);
RHO    = exp(-lambda * abs(tGrid' - tGrid));
C      = chol(RHO, 'lower');

%% 5. Monte Carlo pricing
Nsim    = 100000;
payoffs = priceMC_ExoticCap(nu, B0_T, delta, L0, C, RHO, N, Nsim);

price   = mean(payoffs);
std_err = std(payoffs) / sqrt(Nsim);
fprintf('Exotic Cap price (BMM): %.8f\n', price);
fprintf('Standard error:         %.8f\n', std_err);
fprintf('95%% CI: [%.8f, %.8f]\n', price - 1.96*std_err, price + 1.96*std_err);