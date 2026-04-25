clear; clc; close all;
%% Bootstrap
formatDate = 'dd/mm/yyyy';
[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet);

%% Exercise 1 - Equity Protection Certificate
% data - contract
t0       = datesSet.settlement;       % start date 
P        = 0.95;                      % protection level
alpha    = 1.10;                      % participation coefficient
s        = 0.0130;                    % 130 bp spread over Euribor 3m
N        = 100e6;                     % notional
tN_raw   = addtodate(t0, 5, 'year');   % maturity: 5y after start date, Following Business Day convention
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
fprintf('Deterministic components (per unit of N):\n  floaterPV     = %.6f\n  spreadBPV     = %.6f\n  protectionPV  = %.6f\n\n', floaterPV, spreadBPV, protectionPV);
fprintf('Basket call price (MC, M = %.0e):\n  C_basket = %.6f   (SE = %.6f)\n  95%% CI   = [%.6f , %.6f]\n\n', M, C_basket, SE_basket, CI_basket(1), CI_basket(2));
fprintf('Fair upfront X%%:\n  X       = %.8f %%   (SE = %.6f %%)\n  95%% CI   = [%.6f %% , %.6f %%]\n', 100*X, 100*SE_X, 100*CI_X(1), 100*CI_X(2));
fprintf('  cash    = %.8f M EUR   (on N = %.0f M EUR)\n', N*X/1e6, N/1e6);

%% EXERCISE 2
load('eurostoxx_Poli.mat');   % Load the dataset
strikes = double(cSelect.strikes);   % [31x1] vector of strikes
T = double(cSelect.maturity);   % Time to maturity in years
S0 = double(cSelect.reference);   % ATM Spot price = 3771.06
divYield = double(cSelect.dividends);   % Continuous dividend yield
smiles_mkt = double(cSelect.surface);   % [31x1] implied volatility smile

valuationDate = datenum('15/02/2008', formatDate);
maturityDate = valuationDate + 365;   % Act/365, T=1y

r = interp1(dates, zeroRates, maturityDate, "linear", "extrap");
B = exp(-r * T);         % Discount factor B(0, T) = exp(-rT)
F0 = S0 * exp((r - divYield) * T); %  ATM Forward price
Notional      = 10e6;     %Contract parameters
digitalPayoff = 0.05 * Notional;      % 5% of Notional (500,000 EUR) (paid if S_T > K)
K_digital = F0;      % strike of the digital

% Interpolate smile at the ATM strike K = F0 to obtain ATM implied volatility
sigma_ATM = interp1(strikes, smiles_mkt, K_digital, 'spline');  %we choose spline because we need differentiability

% Black-Scholes price of a cash-or-nothing digital CALL (NO SMILE)
d1_black = ( log(F0 / K_digital) + 0.5 * sigma_ATM^2 * T ) / (sigma_ATM * sqrt(T));
d2_black = d1_black - sigma_ATM * sqrt(T);
% Price of digital = discounted expected payoff under risk-neutral measure
price_Black = B* digitalPayoff * normcdf(d2_black)

%  smile-adjusted pricing; A cash-or-nothing digital can be replicated as the limit of a call spread
eps = 1; 
K_low = K_digital - eps;
K_hi  = K_digital + eps;
sigmas = interp1(strikes, smiles_mkt, [K_low, K_hi], 'spline'); % volatility interpolation

%  Call prices with specific volatility of smile
[C_low, ~] = blkprice(F0, K_low, r, T, sigmas(1));
[C_hi,  ~] = blkprice(F0, K_hi,  r, T, sigmas(2));
price_Smile = digitalPayoff * (C_low - C_hi) / (2 * eps) %  Digital price with Smile (centered differences)

diff_EUR  = price_Smile - price_Black; %Comparison
diff_perc = 100 * diff_EUR / price_Black;
%% PLOT THE IMPLIED VOLATILITY SMILE
% Extended fine grid for smooth smile plot
K_fine = linspace(strikes(1), strikes(end), 500);
vol_fine = interp1(strikes, smiles_mkt, K_fine, 'spline');
figure('Name', 'Implied Volatility Smile - Eurostoxx 50 (15 Feb 2008)');
plot(strikes, smiles_mkt*100, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6); hold on;
plot(K_fine, vol_fine*100, 'b-', 'LineWidth', 1.5);
xline(F0, 'r--', 'LineWidth', 1.5, 'Label', sprintf('ATM F_0 = %.1f', F0));
ylabel('Implied Volatility (%)');
xlabel('Strike K');
title('Eurostoxx 50 Implied Volatility Smile - 1Y, 15 Feb 2008');
legend('Market vols', 'Spline interpolation', 'ATM price F0', 'Location', 'NorthEast');
grid on;
%% EXERCISE 5 - Global Calibration of nMV model
alpha = 2/3;
C_mkt = BS_call(F0, strikes, B, smiles_mkt, T); % Convert market implied vols into market call prices (not IVs) because the calibration minimizes price errors
params0 = [0.20, 2.0, 4.0]; % initial guess  [sigma, kappa, eta]
lb = [1e-4, 1e-4, -50.0];    % lower bounds
ub = [2.00, 10.0,  50.0];    % upper bounds

% Define the constraint eta > -omega_bar as a nonlinear inequality constraint for fmincon
nonlcon = @(p) deal( -p(3) - (1-alpha)/(p(2)*p(1)^2) , []);

% Objective function for global calibration: minimize sum (Cmodel-Cmkt)^2
obj = @(p) sum( (nMV_call_FFT(p, alpha, F0, strikes, B, T) - C_mkt).^2 );

%  Run the optimization with fmincon
options = optimoptions('fmincon', 'Display', 'iter','MaxIterations', 5000, 'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10);
[p_opt, fval] = fmincon(obj, params0, [], [], [], [], lb, ub, nonlcon, options);

sigma_opt = p_opt(1);  kappa_opt = p_opt(2);  eta_opt   = p_opt(3);
fprintf('\n--- Calibrated Parameters (alpha = 2/3) ---\n  sigma = %.6f\n  kappa = %.6f\n  eta   = %.6f\n  Objective value: %.6e\n', ...
    sigma_opt, kappa_opt, eta_opt, fval);   %print results

C_model  = nMV_call_FFT(p_opt, alpha, F0, strikes, B, T);
IV_model = zeros(size(strikes));    %model implied volatilities at market strikes
for i = 1:length(strikes)
    % blkimpv inverts the Black formula: given call price, find sigma
    IV_model(i) = blkimpv(F0, strikes(i), r, T, C_model(i));
end
% RMSE in implied volatility
RMSE = sqrt(mean((IV_model - smiles_mkt).^2)) * 100;
fprintf('RMSE = %.4f%%\n', RMSE);
%% Plot market vs model implied volatility
figure;
x_mkt = log(F0 ./ strikes);
plot(x_mkt, smiles_mkt*100, 'b o-', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', 'Market');
hold on;
plot(x_mkt, IV_model*100, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('nMV \\alpha=2/3 (RMSE=%.3f%%)', RMSE));
xlabel('Moneyness x = ln(F_0/K)');
ylabel('Implied Volatility (%)');
title('Global Calibration – nMV \alpha=2/3');
legend('Location', 'best');
grid on;