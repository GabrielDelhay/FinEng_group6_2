%% 
%  runPricingFourier_ex3.m  
%%

clear; clc; close all;

%% Bootstrap
formatDate = 'dd/mm/yyyy';
today = datenum('15 Feb 08');
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
r = interp1(dates, zeroRates, today);                      
B = exp(-r * T);
F = S0 * exp((r - q) * T);      % ATM forward

fprintf('=== Market Data ===\n');
fprintf('S0=%.4f  F=%.4f  B=%.6f  r=%.4f  q=%.4f  T=%.1f\n\n',...
    S0, F, B, r, q, T);

%% parameters
p_plus  = 1.5;
p_minus = 0.9;

% Martingale condition: phi^c(-i) = 1  =>  E[e^X] = 1
%   phi^c(-i) = exp(mu) / ((1+1/p+)(1-1/p-))  = 1
%   => mu = -log((1+1/p+)(1-1/p-))
mu = -log((1 + 1/p_plus) * (1 - 1/p_minus));
fprintf('Martingale drift: mu = %.8f\n\n', mu);

% %% ---- 3. TARGET MONEYNESS ---------------------------------
% % x = log(F/K),  three targets from the exercise
% x_vec = [-0.05223, 0, 0.15];
% K_vec = F * exp(-x_vec);
% 
% %% ---- 4. PRICING ------------------------------------------
% C_quad = zeros(1,3);
% C_mc   = zeros(1,3);
% C_mc_se= zeros(1,3);
% 
% for j = 1:3
%     x = x_vec(j);
%     K = K_vec(j);
% 
%     % (a) Quadrature
%     I_q      = integralLewis_Quadrature(x, p_plus, p_minus, mu);
%     C_quad(j) = B * (F - sqrt(F*K)/pi * I_q);
% 
%     % (b) Residues  [computed inside same quadrature call - see function]
%     %     -> same result, consistency check printed inside function
% 
% end
% 
% % (c) Monte Carlo
% [C_mc, C_mc_se] = priceMC(x_vec, F, B, p_plus, p_minus, mu, 1e6);
% 
% % (d) FFT  -  two parameter choices
% [C_fft1, params1] = integralLewis_FFT(x_vec, F, B, p_plus, p_minus, mu, ...
%     'M', 15, 'dz', []);   % dz deduced from market strikes
% [C_fft2, params2] = integralLewis_FFT(x_vec, F, B, p_plus, p_minus, mu, ...
%     'M', 12, 'dz', []);   % coarser M, same deduction logic
% 
% %% ---- 5. RESULTS TABLE ------------------------------------
% fprintf('\n%s\n', repmat('=',1,62));
% fprintf('%-16s  %12s  %12s  %12s\n','Method','x=-5.223%','x=0%(ATM)','x=+15%');
% fprintf('%s\n', repmat('-',1,62));
% fprintf('%-16s  %12.6f  %12.6f  %12.6f\n','Quadrature', C_quad);
% fprintf('%-16s  %12.6f  %12.6f  %12.6f\n','MC (1M paths)', C_mc);
% fprintf('%-16s  %12.6f  %12.6f  %12.6f\n','  MC 95%% CI ±', 2*C_mc_se);
% fprintf('%-16s  %12.6f  %12.6f  %12.6f\n', ...
%     sprintf('FFT M=%d',params1.M), C_fft1);
% fprintf('%-16s  %12.6f  %12.6f  %12.6f\n', ...
%     sprintf('FFT M=%d',params2.M), C_fft2);
% fprintf('%s\n\n', repmat('=',1,62));
% 
% %% ---- 6. FFT SENSITIVITY ANALYSIS -------------------------
% fprintf('=== FFT sensitivity: varying M and dz ===\n');
% M_list  = [10, 11, 12, 13, 14, 15];
% dz_list = [0.001, 0.0025, 0.005, 0.01, 0.02];
% 
% % Vary M (dz deduced each time)
% atm_vs_M = zeros(size(M_list));
% for ii = 1:length(M_list)
%     [C_tmp,~] = integralLewis_FFT(0, F, B, p_plus, p_minus, mu, ...
%         'M', M_list(ii), 'dz', []);
%     atm_vs_M(ii) = C_tmp;
% end
% 
% % Vary dz (M=15 fixed)
% atm_vs_dz = zeros(size(dz_list));
% for ii = 1:length(dz_list)
%     [C_tmp,~] = integralLewis_FFT(0, F, B, p_plus, p_minus, mu, ...
%         'M', 15, 'dz', dz_list(ii));
%     atm_vs_dz(ii) = C_tmp;
% end
% 
% figure('Name','FFT Sensitivity','Position',[100 100 1000 420]);
% subplot(1,2,1);
% plot(M_list, atm_vs_M, 'bo-','LineWidth',1.5,'MarkerSize',7);
% yline(C_quad(2),'r--','Quadrature','LineWidth',1.5);
% xlabel('M  (N = 2^M)'); ylabel('ATM Call Price');
% title('Sensitivity to M  (dz deduced)'); grid on;
% 
% subplot(1,2,2);
% semilogx(dz_list, atm_vs_dz, 'bs-','LineWidth',1.5,'MarkerSize',7);
% yline(C_quad(2),'r--','Quadrature','LineWidth',1.5);
% xlabel('dz (moneyness step)'); ylabel('ATM Call Price');
% title('Sensitivity to dz  (M=15 fixed)'); grid on;
% sgtitle('FFT Parameter Sensitivity');
