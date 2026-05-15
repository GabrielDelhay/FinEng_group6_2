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

%% question a) NIG model
[p_opt, alpha_calib, S0, divYield, strikes, smiles_mkt] = ...
    calibrate_NIG(cSelect, dates, zeroRates, formatDate);

% Contract parameters
N_cert = 100e6;  K_str = 3200;  c1 = 0.06;  c2 = 0.02;  spread = 0.0130;

% Payment dates and coupon dates
adj_follow  = @(d) d.*isbusday(d) + busdate(d,1).*~isbusday(d);
T1          = adj_follow(datemnth(t0, 12));
T2          = adj_follow(datemnth(t0, 24));
T1_reset    = busdate(busdate(T1,-1),-1);
T2_reset    = busdate(busdate(T2,-1),-1);
q_dates     = adj_modfollow(datemnth(t0, (3:3:24)'));
assert(T1 == q_dates(4) && T2 == q_dates(8), 'date mismatch');

% Discount factors
B1  = linearRateInterp(dates, discounts, t0, T1);
B2  = linearRateInterp(dates, discounts, t0, T2);
B_q = linearRateInterp(dates, discounts, t0, q_dates);

% Floating leg BPVs (Act/360)
prev_q  = [t0; q_dates(1:end-1)];
delta_q = yearfrac(prev_q, q_dates, 2);
BPV_y1  = sum(delta_q(1:4) .* B_q(1:4));
BPV_y2  = sum(delta_q(5:8) .* B_q(5:8));

% Forward and log-moneyness at reset date
tau1 = yearfrac(t0, T1_reset, 3);
r1   = interp1(dates, zeroRates, T1_reset, 'linear', 'extrap');
F1   = S0 * exp((r1 - divYield) * tau1);
x_K  = log(K_str / F1);

% Risk-neutral probability and upfront — NIG
p1        = digital_put_NIG(p_opt, alpha_calib, x_K, tau1);
% upfront
X_NIG     = price_certificate(p1, N_cert, c1, c2, B1, B2, BPV_y1, BPV_y2, spread);

fprintf('Upfront NIG   = %.4f%%\n', X_NIG   * 100);

%% question b) Black model
sigma_black = interp1(strikes, smiles_mkt, K_str, 'spline');
p1_black    = digital_put_Black(F1, K_str, sigma_black, tau1);
[X_Black]   = price_certificate(p1_black, N_cert, c1, c2, B1, B2, BPV_y1, BPV_y2, spread);


fprintf('Upfront Black = %.4f%%\n', X_Black  * 100);
fprintf('Error         = %.4f bp\n', (X_NIG - X_Black) * 1e4);

%% quetion d) 3y expiry with NIG model - Monte Carlo

% New dates and discount factors
T3       = adj_follow(datemnth(t0, 36));
T3_reset = busdate(busdate(T3,-1),-1);
B3       = linearRateInterp(dates, discounts, t0, T3);

q_dates3 = adj_modfollow(datemnth(t0, (3:3:36)'));
B_q3     = linearRateInterp(dates, discounts, t0, q_dates3);
prev_q3  = [t0; q_dates3(1:end-1)];
delta_q3 = yearfrac(prev_q3, q_dates3, 2);
BPV_y3   = sum(delta_q3(9:12) .* B_q3(9:12));

% Forward prices at reset dates (forward dynamics)
tau2  = yearfrac(t0, T2_reset, 3);
r2    = interp1(dates, zeroRates, T2_reset, 'linear', 'extrap');
F2    = S0 * exp((r2 - divYield) * tau2);
x_K2  = log(K_str / F2);

delta_tau = tau2 - tau1;   % increment period T1 -> T2 (~1 year)

% Simulate NIG increments
N_sim = 1e6;
rng(42);
Z1 = sim_NIG(N_sim, p_opt, alpha_calib, tau1);       % X(T1)
Z2 = sim_NIG(N_sim, p_opt, alpha_calib, delta_tau);  % X(T2)-X(T1), independent

% Spot prices at each observation date
S_T1 = F1 * exp(Z1);           % S(T1_reset)
S_T2 = F2 * exp(Z1 + Z2);      % S(T2_reset)

% Early redemption logic
ER_T1   = S_T1 < K_str;                    % early redemption at T1
ER_T2   = ~ER_T1 & (S_T2 < K_str);        % survive T1, ER at T2
survive = ~ER_T1 & ~ER_T2;                 % reach T3

% Risk-neutral probabilities (Monte Carlo estimates)
p1_mc     = mean(ER_T1);
p_joint   = mean(ER_T2);
p_survive = mean(survive);

% upfront
X_NIG3 = price_certificate3y(p1_mc, p_joint, p_survive, N_cert, c1, c2, ...
                              B1, B2, B3, BPV_y1, BPV_y2, BPV_y3, spread);

fprintf('p1  (ER at T1)       = %.4f%%\n', p1_mc*100);
fprintf('p12 (ER at T2)       = %.4f%%\n', p_joint*100);
fprintf('p3  (reach T3)       = %.4f%%\n', p_survive*100);
fprintf('Upfront X (3y NIG)   = %.4f%%\n', X_NIG3*100);

%% EX 2

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