%% EXERCISE 2
%  Compare the price of a cash-or-nothing digital call option computed:
%    (1) Under the Black model (flat volatility = ATM vol, no smile)
%    (2) Taking into account the implied volatility smile (call-spread replication)
clear; clc; close all;

% Load the dataset (struct 'cSelect' contains all market info)
load('eurostoxx_Poli.mat');

strikes = double(cSelect.strikes);   % [31x1] vector of strikes
T = double(cSelect.maturity);   % Time to maturity in years
S0 = double(cSelect.reference);   % ATM Spot price = 3771.06
divYield = double(cSelect.dividends);   % Continuous dividend yield
smiles_mkt = double(cSelect.surface);   % [31x1] implied volatility smile

formatDate = 'dd/mm/yyyy';
valuationdate = '15/02/2008';
valuationDate = datenum(valuationdate, formatDate);
maturityDate = valuationDate + 365;   % Act/365, T=1y
[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet);

r = interp1(dates, zeroRates, maturityDate, "linear", "extrap");
B = exp(-r * T);         % Discount factor B(0, T) = exp(-rT)
%  ATM Forward price
F0 = S0 * exp((r - divYield) * T);

%Contract parameters
Notional      = 10e6;
digitalPayoff = 0.05 * Notional;          % 5% of Notional (500,000 EUR) (paid if S_T > K)
K_digital = F0;      % strike of the digital

% Interpolate smile at the ATM strike K = F0 to obtain ATM implied volatility
sigma_ATM = interp1(strikes, smiles_mkt, K_digital, 'spline');  %we choose spline because we need differentiability

% The Black-Scholes price of a cash-or-nothing digital CALL is: (NO SMILE)
% Digital_Black = B * cashAmount * N(d2)
d1_black = ( log(F0 / K_digital) + 0.5 * sigma_ATM^2 * T ) / (sigma_ATM * sqrt(T));
d2_black = d1_black - sigma_ATM * sqrt(T);
% Price of digital = discounted expected payoff under risk-neutral measure
price_Black = B* digitalPayoff * normcdf(d2_black)

%  smile-adjusted pricing; A cash-or-nothing digital can be replicated as the limit of a call spread
eps = 1; 
K_low = K_digital - eps;
K_hi  = K_digital + eps;

% volatility interpolation
sigmas = interp1(strikes, smiles_mkt, [K_low, K_hi], 'spline');

%  Call prices with specific volatility of smile
[C_low, ~] = blkprice(F0, K_low, r, T, sigmas(1));
[C_hi,  ~] = blkprice(F0, K_hi,  r, T, sigmas(2));

% Prezzo Digital con Smile (Derivata centrale)
price_Smile = digitalPayoff * (C_low - C_hi) / (2 * eps)

%Comparison
diff_EUR  = price_Smile - price_Black;
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
% Convert market implied vols into market call PRICES.  We need prices (not IVs) because the calibration minimizes price errors
C_mkt = BS_call(F0, strikes, B, smiles_mkt, T);

% Initial guess and bounds for parameters params = [sigma, kappa, eta]
params0 = [0.20, 2.0, 4.0]; % initial guess
lb = [1e-4, 1e-4, -50.0];    % lower bounds
ub = [2.00, 10.0,  50.0];    % upper bounds

% Define the constraint eta > -omega_bar as a nonlinear inequality constraint for fmincon
nonlcon = @(p) deal( -p(3) - (1-alpha)/(p(2)*p(1)^2) , []);

% Objective function for GLOBAL calibration: minimize sum (Cmodel-Cmkt)^2
obj = @(p) sum( (nMV_call_FFT(p, alpha, F0, strikes, B, T) - C_mkt).^2 );

%  Run the optimization with fmincon
options = optimoptions('fmincon', 'Display', 'iter','MaxIterations', 5000, 'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10);

[p_opt, fval] = fmincon(obj, params0, [], [], [], [], lb, ub, nonlcon, options);

sigma_opt = p_opt(1);
kappa_opt = p_opt(2);
eta_opt   = p_opt(3);
fprintf('\n--- Calibrated Parameters (alpha = 2/3) ---\n');
fprintf('  sigma = %.6f\n', sigma_opt);
fprintf('  kappa = %.6f\n', kappa_opt);
fprintf('  eta   = %.6f\n', eta_opt);
fprintf('  Objective value: %.6e\n', fval);

%Compute model implied volatilities at market strikes
C_model  = nMV_call_FFT(p_opt, alpha, F0, strikes, B, T);
IV_model = zeros(size(strikes));
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