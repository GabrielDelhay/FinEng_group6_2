clear; clc; close all;

%% Common market data (15-Feb-2008, Eurostoxx, T = 1y)
[F0, B, r, T] = get_datas();
fprintf('Market data: F0=%.4f  B=%.6f  r=%.4f  T=%.2f\n\n', F0, B, r, T);


%% =============================================================
%  Exercise 3 - Double-exponential characteristic function model
%  =============================================================
fprintf('===== Exercise 3 =====\n');

p_plus  = 1.5;
p_minus = 0.9;
mu_ex3  = log((1 - 1/p_plus) * (1 + 1/p_minus));    % martingale drift

% characteristic function as handle (model-specific factory)
phi_ex3 = charFun_DoubleExp(p_plus, p_minus);

x_vec3 = [-0.05223, 0, 0.15];

% (a) quadrature
C_quad3 = Quadrature(x_vec3, phi_ex3, B, F0);

% (b) residuals (closed-form-by-residues, model specific)
C_res3 = zeros(length(x_vec3), 1);
for j = 1:length(x_vec3)
    I_r       = integralLewis_Residuals(x_vec3(j), p_plus, p_minus, mu_ex3);
    C_res3(j) = B * F0 * (1 - exp(-x_vec3(j)/2) / (2*pi) * I_r);
end

% (c) Monte Carlo
C_mc3 = priceMC(x_vec3, F0, B, p_plus, p_minus, mu_ex3, 1e6);

% (d) FFT
N3  = 2^12;
dx3 = 0.05;
[xg3, I_fft3] = compute_FFT(phi_ex3, N3, dx3);
I_int3 = interp1(xg3, I_fft3, x_vec3, 'spline');
C_fft3 = B * F0 * (1 - exp(-x_vec3/2) ./ (2*pi) .* I_int3);

print_EX_3( x_vec3, C_quad3, C_res3, C_mc3, C_fft3)
%% =============================================================
%  Exercise 4 - Normal Mean-Variance Mixture (NIG, alpha = 1/2)
%  =============================================================
fprintf('===== Exercise 4 =====\n');

sigma = 0.20;  kappa = 1;  eta = 3;  alpha = 1/2;
x_vec4 = (-0.25:0.01:0.25).';

% characteristic function as handle
phi_ex4 = charFun_NMVM(sigma, kappa, eta, T, alpha);

% (a) FFT - three choices of FFT parameters
M  = 12;
N4 = 2^M;

% Choice 1: small moneyness step (dz = 0.0025)
dz1 = 0.0025;          dx1 = 2*pi / (N4*dz1);
[xg41, If41] = compute_FFT(phi_ex4, N4, dx1);
C_fft41      = B * F0 * (1 - exp(-x_vec4/2) ./ (2*pi) .* ...
                         interp1(xg41, If41, x_vec4, 'spline'));

% Choice 2: extremum in frequency domain (u_max = 500, total length 1000)
dx2 = 1000 / N4;
[xg42, If42] = compute_FFT(phi_ex4, N4, dx2);
C_fft42      = B * F0 * (1 - exp(-x_vec4/2) ./ (2*pi) .* ...
                         interp1(xg42, If42, x_vec4, 'spline'));

% Choice 3: dz = 0.01, perfectly aligned with the 1% moneyness grid
dz3 = 0.01;            dx3 = 2*pi / (N4*dz3);
[xg43, If43] = compute_FFT(phi_ex4, N4, dx3);
C_fft43      = B * F0 * (1 - exp(-x_vec4/2) ./ (2*pi) .* ...
                         interp1(xg43, If43, x_vec4, 'spline'));

% (b) Quadrature
C_quad4 = Quadrature(x_vec4, phi_ex4, B, F0);

% (c) Monte Carlo (simulating IG and Gaussian)
C_mc4 = priceMC_ex4(x_vec4, F0, B, sigma, kappa, eta, T, 1e6, alpha);

% comparison plot
fig = comparison_plot_ex4(x_vec4, C_quad4, C_fft41, C_fft42, C_fft43, C_mc4);


% errors
[err1, err2, err3, errMC] = errors_ex4(C_quad4, C_fft41, C_fft42, C_fft43, C_mc4);
