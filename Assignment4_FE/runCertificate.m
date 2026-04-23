%% 
%  runCertificate.m
%%

clear; clc; close all;

%% Bootstrap
formatDate = 'dd/mm/yyyy';
[datesSet, ratesSet] = readExcelDataOS('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet);

%% Exercise 1 - Equity Protection Certificate

%% data - contract
t0       = datesSet.settlementDate;   % start date 
P        = 0.95;                      % protection level
alpha    = 1.10;                      % participation coefficient
s        = 0.0130;                    % 130 bp spread over Euribor 3m
N        = 100e6;                     % notional

% maturity: 5y after start date, Following Business Day convention
tN_raw   = addtodate(t0, 5, 'year');
tN       = busdate(tN_raw - 1, 1);

% quarterly payment dates, Modified Following  
nPay     = 20;                            
rawPay   = datemnth(t0, 3 * (1:nPay)');
payDates = busdate(rawPay - 1, 1);

delta    = yearfrac([t0; payDates(1:end-1)], payDates, 2);  % day-count fractions Act/360 for the floating leg

% discount factors at payment dates
tau_pay  = yearfrac(t0, payDates, 3);
r_pay    = interp1(dates, zeroRates, payDates);
B_pay    = exp(-r_pay .* tau_pay);
B_tN     = B_pay(end);

%% data - equity model
S1_0    = 12.30;    sigma1 = 0.201;    d1 = 0.032;   % ENI
S2_0    = 22.10;    sigma2 = 0.183;    d2 = 0.029;   % AXA
rho     = 0.49;

% maturity ttm (Act/365) and zero rate at tN
T       = yearfrac(t0, tN, 3);
r       = interp1(dates, zeroRates, tN);

%% deterministic components of the swap NPV (per unit of notional)
floaterPV    = 1 - B_tN;                     % Euribor floater 
spreadBPV    = s * sum(delta .* B_pay);      % 130 bp spread leg
protectionPV = (1 - P) * B_tN;               % (1-P) paid at maturity by A

%% basket call price via Monte Carlo
M = 1e6;
rng(42);
[C_basket, SE_basket, CI_basket] = priceBasketCall_MC(S1_0, S2_0, sigma1, sigma2, d1, d2, rho, r, T, M);

%% fair upfront X%   (NPV of the swap = 0)
X      = floaterPV + spreadBPV + protectionPV - alpha * C_basket;
SE_X   = alpha * SE_basket;
CI_X   = X + 1.96 * SE_X * [-1, 1];

%% results
fprintf('Deterministic components (per unit of N):\n');
fprintf('  floaterPV     = %.6f\n', floaterPV);
fprintf('  spreadBPV     = %.6f\n', spreadBPV);
fprintf('  protectionPV  = %.6f\n\n', protectionPV);

fprintf('Basket call price (MC, M = %.0e):\n', M);
fprintf('  C_basket = %.6f   (SE = %.6f)\n', C_basket, SE_basket);
fprintf('  95%% CI   = [%.6f , %.6f]\n\n', CI_basket(1), CI_basket(2));

fprintf('Fair upfront X%%:\n');
fprintf('  X       = %.4f %%   (SE = %.4f %%)\n', 100*X, 100*SE_X);
fprintf('  95%% CI   = [%.4f %% , %.4f %%]\n', 100*CI_X(1), 100*CI_X(2));
fprintf('  cash    = %.4f M EUR   (on N = %.0f M EUR)\n', N*X/1e6, N/1e6);