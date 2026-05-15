% runAssignment6_Group6
%
%

clc; close all; clear all;

formatDate = 'dd/mm/yyyy';
[datesSet, ratesSet] = readExcelDataOS('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet);
t0 = dates(1);

%% EXERCICE 1 : Certificate pricing
fprintf('\n=== Exercise 1 ===\n');
load('eurostoxx_Poli.mat');

% NIG calibration
[p_opt, alpha_calib, S0, divYield, strikes, smiles_mkt] = ...
    calibrate_NIG(cSelect, dates, zeroRates, formatDate);

% Contract parameters
N_cert = 100e6;  K_str = 3200;  c1 = 0.06;  c2 = 0.02;  spread = 0.0130;
adj_follow = @(d) d.*isbusday(d) + busdate(d,1).*~isbusday(d);

% Schedules
[T1, T1_reset, B1, BPV_y1] = get_schedule(t0, 12, dates, discounts, adj_follow);
[T2, T2_reset, B2, BPV_y2] = get_schedule(t0, 24, dates, discounts, adj_follow);
[T3, ~,        B3, BPV_y3] = get_schedule(t0, 36, dates, discounts, adj_follow);

% Forwards
[F1, tau1, x_K]  = get_forward(t0, T1_reset, S0, divYield, dates, zeroRates, K_str);
[F2, tau2, ~  ]  = get_forward(t0, T2_reset, S0, divYield, dates, zeroRates, K_str);

% a) NIG
p1    = digital_put_NIG(p_opt, alpha_calib, x_K, tau1);
X_NIG = price_certificate(p1, N_cert, c1, c2, B1, B2, BPV_y1, BPV_y2, spread);
fprintf('Upfront NIG   = %.4f%%\n', X_NIG * 100);

% b) Black
sigma_black = interp1(strikes, smiles_mkt, K_str, 'spline');
p1_black    = digital_put_Black(F1, K_str, sigma_black, tau1);
X_Black     = price_certificate(p1_black, N_cert, c1, c2, B1, B2, BPV_y1, BPV_y2, spread);
fprintf('Upfront Black = %.4f%%\n', X_Black * 100);
fprintf('Error         = %.4f bp\n', (X_NIG - X_Black) * 1e4);

% d) 3y Monte Carlo
rng(42);
Z1 = sim_NIG(1e6, p_opt, alpha_calib, tau1);
Z2 = sim_NIG(1e6, p_opt, alpha_calib, tau2 - tau1);
ER_T1   = F1*exp(Z1)      < K_str;
ER_T2   = ~ER_T1 & (F2*exp(Z1+Z2) < K_str);
survive = ~ER_T1 & ~ER_T2;
X_NIG3  = price_certificate3y(mean(ER_T1), mean(ER_T2), mean(survive), ...
              N_cert, c1, c2, B1, B2, B3, BPV_y1, BPV_y2, BPV_y3, spread);
fprintf('p1  (ER at T1)     = %.4f%%\n', mean(ER_T1)*100);
fprintf('p12 (ER at T2)     = %.4f%%\n', mean(ER_T2)*100);
fprintf('p3  (reach T3)     = %.4f%%\n', mean(survive)*100);
fprintf('Upfront X (3y NIG) = %.4f%%\n', X_NIG3*100);

%% EXERCISE 2: Bermudian Swaption Pricing via Hull-White.

[sigma, a, K, N, nStepsPerYear, dt, dx, muHat, lMax, x] = getTreeParameters();

%% Calendar (struct)
cal = getCalendar(N, t0, nStepsPerYear);
M   = cal.tree.M;
% Pre-compute initial DFs at each tree node
B0_node = arrayfun(@(d) linearRateInterp(dates, discounts, t0, d), cal.tree.dates);
% Pre-compute deltaX matrices and tree probabilities
[deltaX_top, deltaX_inner, deltaX_bot] = builddeltaX(dx, lMax);
[p_inner, p_bot, p_top] = treeProbabilities(a, dt, lMax);

%% Anchor subset for the swap (drop the non-call year 1, keep maturity)
exerciseYF   = cal.anchor.yf(2:end);                 % years 2..10
exercise_idx = cal.anchor.idx(2:end);

%% 2.A — Bermudan price via tree
V_terminal = zeros(2*lMax+1, 1);
priceSwaption = bermudanRollback(V_terminal, x, p_inner, p_bot, p_top, dt, dx, cal.tree.yf, B0_node, exerciseYF, exercise_idx, a, sigma, K, lMax, 1);

%% 2.B — Tree sanity check: reprice B(0, T_N)
V_terminal = ones(2*lMax+1, 1);
B0_10_tree = bermudanRollback(V_terminal, x, p_inner, p_bot, p_top, dt, dx, cal.tree.yf, B0_node, exerciseYF, exercise_idx, a, sigma, K, lMax, 0);                          
B0_10 = linearRateInterp(dates, discounts, t0, cal.tree.dates(end));

%% 2.C — Bounds
exec = cal.anchor.dates(2:end);                      % T_2..T_10
upperBound = getUpperSwaption(dates, discounts, t0, exec, K, a, sigma);
callDates    = cal.anchor.dates(2:end-1);            % T_2..T_9
paymentDates = cal.anchor.dates;                     % T_1..T_10
lowerBound = getLowerSwaption(dates, discounts, t0, callDates, paymentDates, K, a, sigma);

%% Recap
printPricingResults(priceSwaption, B0_10, B0_10_tree, upperBound, lowerBound);