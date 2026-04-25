clear; clc; close all;

%% Ex 1 - Equity Protection Certificate
% Load and bootstrap the zero rate curve
formatDate = 'dd/mm/yyyy';
[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet);

% Contract parameters
t0 = datesSet.settlement;
N = 100e6; P = 0.95; alpha = 1.10; s = 0.0130;
tN = busdate(addtodate(t0, 5, 'year') - 1, 1);

% Payment dates and floating leg setup
payDates = busdate(datemnth(t0, 3 * (1:20)') - 1, 1);
delta = yearfrac([t0; payDates(1:end-1)], payDates, 2);
tau_pay = yearfrac(t0, payDates, 3);
B_pay = exp(-interp1(dates, zeroRates, payDates) .* tau_pay);
B_tN = B_pay(end);

% Equity model parameters (ENI and AXA)
S1_0 = 12.30; sigma1 = 0.201; d1 = 0.032;
S2_0 = 22.10; sigma2 = 0.183; d2 = 0.029;
rho = 0.49; T_ex1 = yearfrac(t0, tN, 3);
r_ex1 = interp1(dates, zeroRates, tN);

% Deterministic components of the swap NPV
floaterPV = 1 - B_tN;
spreadBPV = s * sum(delta .* B_pay);
protectionPV = (1 - P) * B_tN;

% Basket call pricing via Monte Carlo simulation
M = 1e6; rng(42);
[C_basket, SE_basket, CI_basket] = priceBasketCall_MC(S1_0, S2_0, sigma1, sigma2, d1, d2, rho, r_ex1, T_ex1, M);

% Fair upfront X% calculation and output
X = floaterPV + spreadBPV + protectionPV - alpha * C_basket;
SE_X = alpha * SE_basket;
CI_X = X + 1.96 * SE_X * [-1, 1];
print_certificate_results(floaterPV, spreadBPV, protectionPV, M, C_basket, SE_basket, CI_basket, X, SE_X, CI_X, N);

%% Ex 2 - Pricing Digital option
% Load Eurostoxx dataset and extract basic parameters
load('eurostoxx_Poli.mat');
strikes = double(cSelect.strikes); T = double(cSelect.maturity);
S0 = double(cSelect.reference); divYield = double(cSelect.dividends);
smiles_mkt = double(cSelect.surface);
maturityDate = datenum('15/02/2008', formatDate) + 365;

% Compute Forward ATM and discount factor
r = interp1(dates, zeroRates, maturityDate, "linear", "extrap");
B = exp(-r * T); F0 = S0 * exp((r - divYield) * T);
Notional = 10e6; digitalPayoff = 0.05 * Notional; K_digital = F0;

% Standard Black-Scholes pricing (no smile)
sigma_ATM = interp1(strikes, smiles_mkt, K_digital, 'spline');
d1_black = (log(F0/K_digital) + 0.5*sigma_ATM^2*T) / (sigma_ATM*sqrt(T));
d2_black = d1_black - sigma_ATM*sqrt(T);
price_Black = B * digitalPayoff * normcdf(d2_black)

% Smile-adjusted pricing via Call Spread (centered differences)
eps = 1; K_low = K_digital - eps; K_hi = K_digital + eps;
sigmas = interp1(strikes, smiles_mkt, [K_low, K_hi], 'spline');
[C_low, ~] = blkprice(F0, K_low, r, T, sigmas(1));
[C_hi, ~] = blkprice(F0, K_hi, r, T, sigmas(2));
price_Smile = digitalPayoff * (C_low - C_hi) / (2 * eps)
diff_perc = 100 * (price_Smile - price_Black) / price_Black;
fig1 = plot_smile(strikes, smiles_mkt, F0);

%% Ex 5 - Global Calibration of nMV model
% Convert market IVs to Call prices for calibration objective
alpha_calib = 2/3;
C_mkt = BS_call(F0, strikes, B, smiles_mkt, T);

% Optimization setup: initial guess, bounds, and nonlinear constraints
params0 = [0.20, 2.0, 4.0];
lb = [1e-4, 1e-4, -50.0]; ub = [2.00, 10.0, 50.0];
nonlcon = @(p) deal(-p(3) - (1-alpha_calib)/(p(2)*p(1)^2), []);
obj = @(p) sum((nMV_call_FFT(p, alpha_calib, F0, strikes, B, T) - C_mkt).^2);

% Run fmincon to minimize price errors
options = optimoptions('fmincon', 'Display', 'iter','MaxIterations', 5000, 'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10);
[p_opt, fval] = fmincon(obj, params0, [], [], [], [], lb, ub, nonlcon, options);
fprintf('\n--- Calibrated Params (alpha=2/3) ---\n  sigma = %.6f\n  kappa = %.6f\n  eta   = %.6f\n  Obj: %.6e\n', p_opt(1), p_opt(2), p_opt(3), fval);

% Compute model IVs from calibrated prices and evaluate RMSE
C_model = nMV_call_FFT(p_opt, alpha_calib, F0, strikes, B, T);
IV_model = zeros(size(strikes));
for i = 1:length(strikes)
    IV_model(i) = blkimpv(F0, strikes(i), r, T, C_model(i));
end
RMSE = sqrt(mean((IV_model - smiles_mkt).^2)) * 100;
fprintf('RMSE = %.4f%%\n', RMSE);
fig2 = plot_calibration(F0, strikes, smiles_mkt, IV_model, RMSE);