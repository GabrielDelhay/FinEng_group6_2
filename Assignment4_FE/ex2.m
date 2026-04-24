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
price_Black = B * digitalPayoff * normcdf(d2_black)

%  SMILE-ADJUSTED PRICING (call-spread replication)
%
%  A cash-or-nothing digital can be replicated as the limit of a call spread:
%   Digital(K) = lim_{eps -> 0}  [ C(K - eps) - C(K + eps) ] / (2eps)
%
%  In discrete form with finite epsilon:
%    Digital(K) ≈ [ C(K - eps) - C(K + eps) ] / (2eps)
eps = 0.1; 
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
xline(K_digital, 'k:', 'LineWidth', 1.2);
ylabel('Implied Volatility (%)');
xlabel('Strike K');
title('Eurostoxx 50 Implied Volatility Smile - 1Y, 15 Feb 2008');
legend('Market vols', 'Spline interpolation', 'ATM Forward', 'Location', 'NorthEast');
grid on;



%% EXERCISE 5 - Global Calibration of nMV model

alpha = 2/3; 
% Convert market implied vols into market call PRICES.  We need prices (not IVs) because the calibration minimizes price errors
C_mkt = BS_call(F0, strikes, B, smiles_mkt, T);

% Define initial guess and bounds for parameters params = [sigma, kappa, eta]
% From slide 20 (Parsimony):
%   sigma : average volatility (must be > 0)
%   kappa : vol-of-vol, controls smile curvature (must be > 0)
%   eta   : skew parameter, controls smile asymmetry. eta=0 => symmetric smile: sigma_B(x) = sigma_B(-x)
params0 = [0.20, 2.0, 4.0]; % initial guess
lb = [1e-4, 1e-4, 0.0];    % lower bounds
ub = [2.00, 10.0,  50.0];    % upper bounds


% Define the constraint eta > -omega_bar
%   The characteristic function phi(xi) must be analytic in the strip
%   that contains both the normalization point (Im=0) and the martingale property point (Im=-1).
%   This requires: eta > -omega_bar
%   where omega_bar = (1-alpha)/(kappa*sigma^2)  
% We write this as a nonlinear inequality constraint for fmincon
nonlcon = @(p) deal( -p(3) - (1-alpha)/(p(2)*p(1)^2) , []);

% Define the objective function for GLOBAL calibration
%   d = sum_i  omega_i * |c(x_i, t_i; p) - c_i|^2obj = @(p) sum( (nMV_call_FFT(p, alpha, F0, strikes, B, T) - C_mkt).^2 );

obj = @(p) sum( (nMV_call_FFT(p, alpha, F0, strikes, B, T) - C_mkt).^2 );

%  Run the optimization with fmincon
options = optimoptions('fmincon', ...
    'Display',             'iter', ...
    'MaxIterations',       5000,   ...
    'OptimalityTolerance', 1e-10,  ...
    'StepTolerance',       1e-10);

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
% We have model call prices from FFT, now we invert the Black formula
% to get implied vols and compare with market
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
plot(x_mkt, IV_model*100, 'r-', 'LineWidth', 2, ...
    'DisplayName', sprintf('nMV \\alpha=2/3 (RMSE=%.3f%%)', RMSE));
xlabel('Moneyness x = ln(F_0/K)');
ylabel('Implied Volatility (%)');
title('Global Calibration – nMV \alpha=2/3, Eurostoxx50 15-Feb-2008');
legend('Location', 'best');
grid on;
%%

function lnL = laplace_exponent(omega, p, alpha, T)
    % Laplace exponent of the tempered stable subordinator G (slide 15)
    % ln L[omega] = (T/kappa) * (1-alpha)/alpha * [1 - (1 + omega*kappa*sigma^2/(1-alpha))^alpha]
    %
    % Inputs:
    %   omega : can be a scalar or vector (real or complex)
    %   p     : [sigma, kappa, eta]
    %   alpha : stability index (2/3 in our case)
    %   T     : time to maturity (= Delta_t in slide notation)
    
    sigma = p(1);
    kappa = p(2);
    
    lnL = (T / kappa) .* (1 - alpha) / alpha .* ...
          (1 - (1 + omega .* kappa .* sigma^2 ./ (1 - alpha)).^alpha);
end

function phi = char_fun(xi, p, alpha, T)
    % Characteristic function of the log-return ft = ln(Ft/F0) (slide 15)
    %
    % phi(xi) = exp{ -i*xi * ln L[eta] } * L[ (xi^2 + i*(1+2*eta)*xi) / 2 ]
    %
    % where:
    %   - exp{-i*xi * ln L[eta]} is the NORMALIZATION term (slide 17)
    %     it ensures E[e^{ft}] = 1 (martingale property)
    %   - L[(xi^2 + i*(1+2*eta)*xi)/2] is the main Laplace transform
    %     evaluated at a COMPLEX argument (this is why we need analyticity)
    %
    % Note: ln L[eta] is real since eta is real
    
    eta = p(3);
    
    % Normalization: ln L evaluated at real point eta (slide 15, 17)
    lnL_eta = laplace_exponent(eta, p, alpha, T);
    
    % Main argument: complex number (xi^2 + i*(1+2*eta)*xi) / 2  (slide 15)
    omega_complex = (xi.^2 + 1i .* (1 + 2*eta) .* xi) / 2;
    
    % Laplace exponent at complex argument
    lnL_main = laplace_exponent(omega_complex, p, alpha, T);
    
    % Final characteristic function (slide 15)
    phi = exp(-1i .* xi .* lnL_eta + lnL_main);
end

function C = nMV_call_FFT(p, alpha, F0, K_vec, B, T)
    % Lewis (2001) via FFT  –  slide 10, 12, 13
    %
    % Integral:  I(x) = int_{-inf}^{+inf} dxi/(2pi) * e^{-i*xi*x}
    %                   * phi(-xi-i/2) / (xi^2 + 1/4)
    %
    % Discretization (slide 13):
    %   xi_j = xi_1 + (j-1)*dxi,   j = 1,...,N
    %   x_k  = x_1  + (k-1)*dx
    %   dxi * dx = 2*pi/N
    
    M   = 14;
    N   = 2^M;          % 4096 points
    dxi = 0.01;         % frequency step  (tune this)
    dx  = 2*pi/(N*dxi); % moneyness step  (forced by constraint)
    
    % xi grid: centered around 0  (xi goes from negative to positive)
    j    = (0:N-1)';
    xi1  = -(N/2)*dxi;           % starting frequency
    xi_j = xi1 + j*dxi;         % xi_j = xi_1 + (j-1)*dxi
    
    % x grid: centered around 0
    x1   = -(N/2)*dx;
    x_k  = x1 + j*dx;
    
    % Integrand evaluated on xi grid (this is what we FFT-transform)
    % phi must be analytic in strip Im(xi) in [-1,0]: OK by construction
    phi_vals  = char_fun(-xi_j - 1i/2, p, alpha, T);   % phi(-xi-i/2)
    integrand = phi_vals ./ (xi_j.^2 + 0.25);           % divide by xi^2+1/4
    
    % Slide 13: to use FFT we write
    %   I(x_k) = (dxi/2pi) * sum_j integrand(xi_j) * e^{-i*xi_j*x_k}
    %          = (dxi/2pi) * e^{-i*xi_1*x_k} * FFT[ e^{-i*(j-1)*dxi*x_1} * integrand_j ]
    %
    % Define f_j (the sequence fed into FFT):
    f_j = integrand .* exp(-1i * xi_j * x1);  % multiply by e^{-i*xi_j*x_1}
    % Note: xi_j*x1 = (xi_1 + (j-1)*dxi)*x1
    %               = xi_1*x1  +  (j-1)*dxi*x1
    % so f_j = integrand * e^{-i*xi_1*x1} * e^{-i*(j-1)*dxi*x1}
    % the e^{-i*xi_1*x1} is a global phase (const), absorbed in prefactor
    
    % FFT
    FFT_out = fft(f_j);   % FFT_out(k) = sum_j f_j * e^{-2pi*i*(j-1)*(k-1)/N}
    
    % Prefactor to convert FFT output to the integral at x_k  (slide 13)
    % I(x_k) = (dxi/2pi) * e^{-i*xi_1*x_k} * FFT_out(k)
    prefactor = exp(-1i * xi1 * x_k);
    I_k = (dxi / (2*pi)) * prefactor .* FFT_out;
    
    % Lewis formula: c(x)/(B*F0) = 1 - e^{-x/2} * Re{ I(x) }
    C_grid = B * F0 * (1 - exp(-x_k/2) .* real(I_k));
    C_grid = max(C_grid, 0);
    
    % Interpolate at requested strikes
    x_req = log(F0 ./ K_vec(:));
    
    % Check all requested moneyness values fall inside the grid
    if any(x_req < x_k(1)) || any(x_req > x_k(end))
        warning('Some strikes fall outside FFT moneyness grid! Increase M or adjust dxi.');
    end
    
    C = interp1(x_k, C_grid, x_req, 'spline');
    C = max(real(C), 0);
end

function C = BS_call(F0, K, B, sigma, T)
    % Black (1976) call price formula (slide 6)
    % c(K) = B(t0,t) * { F0*N[d1] - K*N[d2] }
    % with d_{1,2} = ln(F0/K)/sqrt(T*sigma^2) +/- (1/2)*sqrt(T*sigma^2)
    d1 = (log(F0./K) + 0.5*sigma.^2.*T) ./ (sigma.*sqrt(T));
    d2 = d1 - sigma.*sqrt(T);
    C  = B .* (F0 .* normcdf(d1) - K .* normcdf(d2));
end
