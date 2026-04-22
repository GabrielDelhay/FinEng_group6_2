%% 
%  runPricingFourier_ex3.m  
%%

clear; clc; close all;

%% Bootstrap
formatDate = 'dd/mm/yyyy';
maturity = datenum('15-Feb-2009');
[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet);

%% Exercise 4
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
x_vec = [-0.05223, 0, 0.15];     % x = log(F0/K)

sigma = 0.20;  kappa = 1;  eta = 3;  alpha = 1/2;
x_vec = (-0.25:0.01:0.25).';

for j = 1:length(x_vec)
    x = x_vec(j);% Quadrature method
    Integrand = @(u) real(lewisIntegrand_ex4(u, x, sigma, kappa, eta, T, alpha));
    I_q = 2 * quadgk(Integrand, 0, Inf); % Integrand is even in the complex sense
    C_quad(j) = B * F0 * (1 - exp(-x/2) / (2*pi) * I_q);
end
% Compute FFT
N = 2^12; du = 0.05;
[x_grid, I_fft]  = compute_FFT_ex4(sigma, kappa, eta, T, alpha, N, du);
x_targets = x_vec;                      % la tua griglia -25% .. +25%
I_interp = interp1(x_grid, I_fft, x_targets, 'spline');
C_fft = B * F0 * (1 - exp(-x_targets/2) ./ (2*pi) .* I_interp);
% Compute MC
C_mc_ = priceMC_ex4(x_vec, F0, B, sigma, kappa, eta, T, 1e6, alpha);


%% results

fprintf('Quadrature:  %.4f\n', C_quad);
fprintf('MonteCarlo:  %.4f\n', C_mc_);
fprintf('FFT:  %.4f\n', C_fft);
