% runAssignment6_Group6
%
%
clc; close all; clear all;

formatDate = 'dd/mm/yyyy';
[datesSet, ratesSet] = readExcelDataOS('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet);
t0 = dates(1);  % Settlement date: 19/02/2008



%% EXERCISE 1 : Certificate Pricing 

% data
load('eurostoxx_Poli.mat');

%% Global Calibration of NIG model

strikes = double(cSelect.strikes); 
T = double(cSelect.maturity);
S0 = double(cSelect.reference); 
divYield = double(cSelect.dividends);
smiles_mkt = double(cSelect.surface);

maturityDate = datenum('15/02/2008', formatDate) + 365;
r = interp1(dates, zeroRates, maturityDate, "linear", "extrap");
B = exp(-r * T);
F0 = S0 * exp((r - divYield) * T);

% Convert market IVs to Call prices for calibration objective
alpha_calib = 1/2; %NIG
C_mkt = BS_call(F0, strikes, B, smiles_mkt, T);

% Optimization setup: initial guess, bounds, and nonlinear constraints
params0 = [0.20, 2.0, 4.0];
lb = [1e-4, 1e-4, -50.0]; ub = [2.00, 10.0, 50.0];
nonlcon = @(p) deal(-p(3) - (1-alpha_calib)/(p(2)*p(1)^2), []);
obj = @(p) sum((nMV_call_FFT(p, alpha_calib, F0, strikes, B, T) - C_mkt).^2);

% Run fmincon to minimize price errors
options = optimoptions('fmincon', 'Display', 'iter','MaxIterations', 5000, 'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10);
[p_opt, fval] = fmincon(obj, params0, [], [], [], [], lb, ub, nonlcon, options);
fprintf('\n--- Calibrated Params (alpha=1/2) ---\n  sigma = %.6f\n  kappa = %.6f\n  eta   = %.6f\n  Obj: %.6e\n', p_opt(1), p_opt(2), p_opt(3), fval);

% Compute model IVs from calibrated prices and evaluate RMSE
C_model = nMV_call_FFT(p_opt, alpha_calib, F0, strikes, B, T);
IV_model = zeros(size(strikes));
for i = 1:length(strikes)
    IV_model(i) = blkimpv(F0, strikes(i), r, T, C_model(i));
end
RMSE = sqrt(mean((IV_model - smiles_mkt).^2)) * 100;
fprintf('RMSE = %.4f%%\n', RMSE);

% certificate 
N_cert = 100e6;    % 100 MIO EUR notional
spread = 0.0130;   % Party A pays Euribor 3m + 1.30% (Act/360, quarterly)
K_str  = 3200;     % STOXX50 strike
c1     = 0.06;     % Year 1 coupon (if STOXX < K)
c2     = 0.02;     % Year 2 coupon (fixed, last year)

% Payment dates

% Following BDC: keep if already a business day, else next BD
adj_follow = @(d) d .* isbusday(d) + busdate(d, 1) .* ~isbusday(d);

% Annual coupon dates — Following BDC, 30/360 (bond & swap coupon leg)
T1 = adj_follow(datemnth(t0, 12));   % ~19 Feb 2009
T2 = adj_follow(datemnth(t0, 24));   % ~19 Feb 2010

% Coupon Reset Dates = 2 BD before respective coupon payment date (in arrears)
% These are the STOXX50 observation dates used in the NIG model
T1_reset = busdate(busdate(T1, -1), -1);   % 2 BD before T1
T2_reset = busdate(busdate(T2, -1), -1);   % 2 BD before T2

% Quarterly dates — Modified Following, Act/360 (swap floating leg)
q_dates = adj_modfollow(datemnth(t0, (3:3:24)'));   % 8 quarterly dates

% T1/T2 coincide with Q4/Q8 (annual coupons on same schedule as quarterly)
assert(T1 == q_dates(4) && T2 == q_dates(8), 'Quarterly/annual date mismatch');

% Discount factors from bootstrap
B1 = linearRateInterp(dates, discounts, t0, T1);
B2 = linearRateInterp(dates, discounts, t0, T2);
B_q = linearRateInterp(dates, discounts, t0, q_dates);

% Floating leg BPVs (Act/360 daycount)
prev_q  = [t0; q_dates(1:end-1)];
delta_q = yearfrac(prev_q, q_dates, 2);    % Act/360
BPV_y1  = sum(delta_q(1:4) .* B_q(1:4));  % year 1 quarters
BPV_y2  = sum(delta_q(5:8) .* B_q(5:8));  % year 2 quarters

% Floating leg NPV (standard floater decomposition)
% NPV(Euribor flat) = N*(B_start - B_end); spread = N*s*BPV
NPV_float_y1 = N_cert * (1 - B1) + spread * N_cert * BPV_y1;
NPV_float_y2 = N_cert * (B1 - B2) + spread * N_cert * BPV_y2;

% Risk-neutral probability p1 = Q(S(T1_reset) < K) under NIG
% tau1 uses the Reset Date (observation of STOXX50), not the payment date
tau1 = yearfrac(t0, T1_reset, 3);   % Act/365
r1   = interp1(dates, zeroRates, T1_reset, 'linear', 'extrap');
F1   = S0 * exp((r1 - divYield) * tau1);   % forward at T1_reset

% x_K = log(K/F1): log-moneyness for the binary put observation
x_K  = log(K_str / F1);

% Gil-Pelaez CDF inversion:  Q(f < x_K) = 1/2 + 1/pi * int_0^inf Im[e^{-iux} phi(u)] / u du
phi_T1       = @(u) char_fun(u, p_opt, alpha_calib, tau1);
gp_integrand = @(u) imag(exp(-1i * u * x_K) .* phi_T1(u)) ./ u;
p1 = 0.5 - (1/pi) * integral(gp_integrand, 1e-8, 500, 'AbsTol', 1e-10, 'RelTol', 1e-8);

% NPV of equity coupons (Party B pays annually)
% At T1: 6%*N if STOXX < K  (early redemption case)
% At T2: 2%*N only if no early redemption (STOXX >= K at T1)
NPV_coupons = N_cert * (c1 * B1 * p1 + c2 * B2 * (1 - p1));

% Total floating NPV with early redemption
% Year 1 always paid; year 2 only if no ER (prob 1-p1)
NPV_float = NPV_float_y1 + (1 - p1) * NPV_float_y2;

% Upfront X% (mid-market: NPV_swap = 0 for Party B)
% Party B: pays X*N + equity coupons, receives Euribor+spread
% => X*N = NPV_float - NPV_coupons
X_upfront = (NPV_float - NPV_coupons) / N_cert;

fprintf('S0 = %.2f  |  F1 = %.2f  |  K = %.0f\n', S0, F1, K_str);
fprintf('x_K = log(K/F1) = %.6f\n', x_K);
fprintf('p1 = Q(STOXX50 < 3200) = %.4f%%\n', p1 * 100);
fprintf('B1 = %.6f  |  B2 = %.6f\n', B1, B2);
fprintf('NPV float (with ER)  = %14.2f EUR\n', NPV_float);
fprintf('NPV equity coupons   = %14.2f EUR\n', NPV_coupons);
fprintf('Upfront X            =     %.4f%%\n', X_upfront * 100);
fprintf('Upfront X * N        = %14.2f EUR\n', X_upfront * N_cert);

% question b) Black model 
% Implied vol at K = 3200 (interpolate on the market smile)
sigma_black = interp1(strikes, smiles_mkt, K_str, 'spline');

% Black digital put: Q(S(T1) < K) = N(-d2)
d2_black = (log(F1 / K_str) - 0.5 * sigma_black^2 * tau1) / (sigma_black * sqrt(tau1));
p1_black = normcdf(-d2_black);

% Same NPV structure, just replace p1 with p1_black
NPV_coupons_black = N_cert * (c1 * B1 * p1_black + c2 * B2 * (1 - p1_black));
NPV_float_black   = NPV_float_y1 + (1 - p1_black) * NPV_float_y2;
X_upfront_black   = (NPV_float_black - NPV_coupons_black) / N_cert;

fprintf('sigma_black (at K=3200) = %.4f%%\n', sigma_black * 100);
fprintf('p1_black = Q(STOXX50 < 3200) = %.4f%%\n', p1_black * 100);
fprintf('Upfront X (Black)       =     %.4f%%\n', X_upfront_black * 100);
fprintf('Delta NIG vs Black      =     %.4f bp\n', (X_upfront - X_upfront_black) * 1e4);

% Calculate the difference in upfront costs between NIG and Black models
delta_upfront = (X_upfront - X_upfront_black) * 1e4;  % in basis points
fprintf('Delta Upfront (NIG - Black) = %.4f bp\n', delta_upfront);