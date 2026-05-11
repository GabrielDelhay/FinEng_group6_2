% runAssignment6_Group6
%
%
clc; close all; clear all;

formatDate = 'dd/mm/yyyy';
[datesSet, ratesSet] = readExcelDataOS('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet);
t0 = dates(1);  % Settlement date: 19/02/2008


%% EXERCISE 1

%% Ex 5 - Global Calibration of nMV model
load('eurostoxx_Poli.mat');
strikes = double(cSelect.strikes); 
T = double(cSelect.maturity);
S0 = double(cSelect.reference); 
divYield = double(cSelect.dividends);
smiles_mkt = double(cSelect.surface);
maturityDate = datenum('15/02/2008', formatDate) + 365;
r = interp1(dates, zeroRates, maturityDate, "linear", "extrap");
B = exp(-r * T);
F0 = S0 * exp((r - divYield) * T);
Notional = 10e6;

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
