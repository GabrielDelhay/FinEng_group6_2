%% 
%  runPricingFourier_ex3.m  
%%

clear; clc; close all;

%% Bootstrap
formatDate = 'dd/mm/yyyy';
maturity = datenum('15-Feb-2009');
[datesSet, ratesSet] = readExcelDataOS('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet);

%% Exercise 3 - Pricing

%% data
load('eurostoxx_Poli.mat');

S0       = cSelect.reference;          % spot
T        = cSelect.maturity;           % 1 year
q        = cSelect.dividends;          % continuous dividend yield
strikes  = double(cSelect.strikes);    % 31 strikes
IV_mkt   = cSelect.surface;           % implied vols (Black-Scholes)

% Discount factor from bootstrap
r = interp1(dates, zeroRates, maturity);                      
B = exp(-r * T);
F0 = S0 * exp((r - q) * T);      % ATM forward
x_vec = [-0.05223, 0, 0.15];     % x = log(F/K)

fprintf('S0=%.4f  F0=%.4f  B=%.6f  r=%.4f  q=%.4f  T=%.1f\n\n',S0, F0, B, r, q, T);

%% parameters
p_plus  = 1.5;
p_minus = 0.9;

% Martingale condition: phi^c(-i) = 1
mu = log((1 - 1/p_plus) * (1 + 1/p_minus));


%% pricing
C_quad = zeros(1,3);
C_residuals = zeros(1,3);

for j = 1:3
    x = x_vec(j);

    % Quadrature method
    Integrand = @(u) real(lewisIntegrand(u,x,p_plus,p_minus,mu));
    I_q      = quadgk(Integrand, -Inf, Inf);
    C_quad(j) = B * F0 * (1 - exp(-x/2) / (2*pi) * I_q);
    % Residuals method
    I_r = integralLewis_Residuals(x,p_plus,p_minus,mu);
    C_residuals(j) = B * F0 * (1 - exp(-x/2) / (2*pi) * I_r);
    % Monte Carlo
    C_mc = priceMC(x_vec, F0, B, p_plus, p_minus, mu, 1e6);
    % FFT
    [x_grid, C_fft] = compute_FFT(p_plus, p_minus, mu);
    x_targets = [-0.05223, 0, 0.15];
    C_interp = interp1(x_grid, C_fft, x_targets, 'spline');
    C_fft = B * F0 * (1 - exp(-x_targets/2) ./ (2*pi) .* C_interp);
end

%% results
fprintf('Quadrature:  %.4f  %.4f  %.4f\n', C_quad);
fprintf('Residuals:  %.4f  %.4f  %.4f\n', C_residuals);
fprintf('MonteCarlo:  %.4f  %.4f  %.4f\n', C_mc);
fprintf('FFT:  %.4f  %.4f  %.4f\n', C_fft);
