% utilities_ex5.m

function val = logLap(w, sigma, k, dt, alpha)

val = dt/k * (1-alpha)/alpha * (1 - (1 + w*k*sigma^2/(1-alpha)).^alpha);

end

function val = phi_NVM(u, params, alpha)

sigma = params(1);
k     = params(2);
dt    = params(3);   % = Delta t
eta   = params(4);

% ln L[eta]  (scalaire réel)
lnL_eta = logLap(eta, sigma, k, dt, alpha);

% ln L[(u^2 + i*(1+2*eta)*u) / 2]  (vecteur complexe)
w       = (u.^2 + 1i*(1 + 2*eta).*u) / 2;
lnL_w   = logLap(w, sigma, k, dt, alpha);

% phi(u) = exp(-i*u * lnL[eta]) * exp(lnL[w])
%        = exp(-i*u * lnL[eta] + lnL[w])
val = exp(-1i.*u * lnL_eta + lnL_w);

end

function C = priceCall_Lewis(x, F0, B, params, alpha)
%% priceCall_Lewis  -  Call price via Lewis formula
%
%  c(x) / (B*F0) = 1 - exp(-x/2)/(2*pi) * integral_{-inf}^{+inf}
%                  exp(-i*xi*x) * phi(-xi - i/2) / (xi^2 + 1/4) dxi
%
%  Inputs:
%    x      - moneyness = log(F0/K)  [scalar]
%    F0     - ATM forward
%    B      - discount factor
%    params - [sigma, k, dt, eta]
%    alpha  - fixed at 2/3
%
%  Output:
%    C      - call price

K = F0 * exp(-x);

% Integrand: phi evaluated at (-xi - i/2) with variable xi
% Using course convention: integrand = exp(-i*xi*x) * phi(-xi-i/2) / (xi^2+1/4)
integrand = @(xi) real( exp(-1i*xi*x) .* phi_NVM(-xi - 1i/2, params, alpha) ...
                        ./ (xi.^2 + 0.25) );

% Integrate over R (course formula with 1/2pi factor)
I = quadgk(integrand, -Inf, Inf, 'RelTol', 1e-8, 'AbsTol', 1e-10);

% Lewis formula
C = B * F0 * (1 - exp(-x/2) / (2*pi) * I);

end

function sigma_iv = impliedVol(C_mkt, F0, K, B, T)
%% impliedVol  -  BS implied volatility via bisection (fzero)
%
%  Finds sigma such that BS_call(sigma) = C_mkt
%
%  Inputs:
%    C_mkt - market/model call price
%    F0    - forward
%    K     - strike
%    B     - discount factor
%    T     - maturity
%
%  Output:
%    sigma_iv - implied volatility

% Lower bound: intrinsic value check
C_intrinsic = max(B*(F0 - K), 0);
if C_mkt <= C_intrinsic + 1e-10
    sigma_iv = 1e-6;
    return;
end

obj = @(s) BS_call(F0, K, B, s, T) - C_mkt;

try
    sigma_iv = fzero(obj, [1e-4, 5], optimset('TolX', 1e-10, 'Display', 'off'));
catch
    sigma_iv = NaN;
end

end

%% ---- Black-Scholes call (forward formulation) ------------
function C = BS_call(F0, K, B, sigma, T)

d1 = (log(F0/K) + 0.5*sigma^2*T) / (sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);
C  = B * (F0*normcdf(d1) - K*normcdf(d2));

end

function err = objectiveCalib(params, x_vec, IV_mkt, F0, B, T, alpha)

sigma = params(1);
k     = params(2);
dt    = params(3);
eta   = params(4);

% Vérification : argument de logLap doit rester dans le domaine
% 1 + eta*k*sigma^2/(1-alpha) > 0
if 1 + eta*k*sigma^2/(1-alpha) <= 0
    err = 1e10;
    return;
end

err = 0;
for j = 1:length(x_vec)
    x = x_vec(j);
    K = F0 * exp(-x);

    C_mod = priceCall_Lewis(x, F0, B, params, alpha);

    if ~isfinite(C_mod) || C_mod <= 0
        err = err + 1;
        continue;
    end

    IV_mod = impliedVol(C_mod, F0, K, B, T);

    if ~isfinite(IV_mod)
        err = err + 1;
        continue;
    end

    err = err + (IV_mod - IV_mkt(j))^2;
end

end







%% =========================================================
%  ex5.m
%%

clear; clc; close all;

%% data
load('eurostoxx_Poli.mat');

S0      = cSelect.reference;
T       = 365/365;
q       = cSelect.dividends;
strikes = double(cSelect.strikes);
IV_mkt  = cSelect.surface;

% Bootstrap
formatDate = 'dd/mm/yyyy';
[datesSet, ratesSet] = readExcelDataOS('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet);
maturity = datenum('15-Feb-2009');
B        = interp1(dates, discounts, maturity);
r        = -log(B) / T;
F0       = S0 * exp((r - q) * T);

fprintf('S0=%.4f  F0=%.4f  B=%.6f  r=%.4f%%  q=%.4f%%\n\n', ...
    S0, F0, B, r*100, q*100);

% Moneyness x = log(F0/K)
x_vec = log(F0 ./ strikes);

%% parameters
alpha   = 2/3;
params0 = [0.20, 1.0, 1.0, 0.05];

sigma0 = params0(1);
k0     = params0(2);

% Contrainte : eta > -omega_bar
omega_bar = (1-alpha) / (k0 * sigma0^2);
fprintf('omega_bar = %.4f  (doit etre > 0)\n', omega_bar);
fprintf('Contrainte eta > %.4f\n', -omega_bar);

lb = [0.01, 0.01, T,   -omega_bar + 1e-4];
ub = [2.00, 10.0, T,    0.49];

%% calibration (global)
opts = optimoptions('fmincon', ...
    'Display',   'iter', ...
    'TolFun',    1e-8,   ...
    'TolX',      1e-8,   ...
    'MaxFunEvals', 5000);

obj = @(p) objectiveCalib(p, x_vec, IV_mkt, F0, B, T, alpha);

params_calib = fmincon(obj, params0, [], [], [], [], lb, ub, [], opts);

fprintf('\n=== Calibrated parameters ===\n');
fprintf('sigma = %.6f\n', params_calib(1));
fprintf('k     = %.6f\n', params_calib(2));
fprintf('dt    = %.6f\n', params_calib(3));
fprintf('eta   = %.6f\n', params_calib(4));

%% COMPUTE MODEL IV
IV_mod = zeros(size(x_vec));
for j = 1:length(x_vec)
    x = x_vec(j);
    K = F0 * exp(-x);
    C = priceCall_Lewis(x, F0, B, params_calib, alpha);
    IV_mod(j) = impliedVol(C, F0, K, B, T);
end

%% RMSE 
RMSE = sqrt(mean((IV_mod - IV_mkt).^2)) * 100;
fprintf('\nRMSE = %.4f%%\n', RMSE);

%% plot
figure('Name', 'NVM Calibration', 'Position', [100 100 800 500]);
plot(strikes, IV_mkt*100, 'ko-', 'LineWidth', 1.5, 'DisplayName', 'Market');
hold on;
plot(strikes, IV_mod*100, 'r-',  'LineWidth', 2,   'DisplayName', sprintf('NVM \\alpha=2/3  (RMSE=%.3f%%)', RMSE));
xlabel('Strike');
ylabel('Implied Volatility (%)');
title('NVM \alpha=2/3 Calibration — Eurostoxx50, 15-Feb-2008');
legend('Location', 'northeast');
grid on;