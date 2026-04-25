clear; clc; close all;

%% Common market data (15-Feb-2008, Eurostoxx, T = 1y)
[F0, B, r, T] = get_datas();
fprintf('Market data: F0=%.4f  B=%.6f  r=%.4f  T=%.2f\n\n', F0, B, r, T);

%% Ex 3 - Double-exponential characteristic function model
fprintf('===== Exercise 3 =====\n');
p_plus = 1.5; p_minus = 0.9;
mu_ex3 = log((1 - 1/p_plus) * (1 + 1/p_minus)); % martingale drift

phi_ex3 = charFun_DoubleExp(p_plus, p_minus);
x_vec3 = [-0.05223, 0, 0.15];

% (a) Quadrature and (b) Residuals
C_quad3 = Quadrature(x_vec3, phi_ex3, B, F0);
C_res3 = zeros(length(x_vec3), 1);
for j = 1:length(x_vec3)
    I_r = integralLewis_Residuals(x_vec3(j), p_plus, p_minus, mu_ex3);
    C_res3(j) = B * F0 * (1 - exp(-x_vec3(j)/2) / (2*pi) * I_r);
end

% (c) Monte Carlo
C_mc3 = priceMC(x_vec3, F0, B, p_plus, p_minus, mu_ex3, 1e6);

% (d) FFT across different configurations
phi_ex3 = @(v) exp(1i*mu_ex3*v) ./ ((1 - 1i*v/p_plus).*(1 + 1i*v/p_minus));
configs = {
    'N=2^12, du=0.05 (ref)', 2^12, 0.05;  'N=2^12, du=0.25',  2^12, 0.25;
    'N=2^12, du=0.01',       2^12, 0.01;  'N=2^8,  du=0.05',  2^8,  0.05;
    'N=2^14, du=0.05',       2^14, 0.05;  'N=2^12, dz=0.005', 2^12, 2*pi/(2^12*0.005);
    'N=2^12, dz=0.10',       2^12, 2*pi/(2^12*0.10);
};

fprintf('%-30s  %10s  %10s  %10s\n', 'Configuration', 'C(x1)', 'C(x2)', 'C(x3)');
fprintf('%s\n', repmat('-',1,65));
for k = 1:size(configs,1)
    Nk = configs{k,2}; duk = configs{k,3};  
    [xg, I_fft] = compute_FFT(phi_ex3, Nk, duk);
    
    if any(x_vec3 < xg(1)) || any(x_vec3 > xg(end))
        warning('Strike off the grid: %s', configs{k,1}); C_fft3 = [NaN NaN NaN];
    else
        I_int = interp1(xg, I_fft, x_vec3, 'spline');
        C_fft3 = B * F0 * (1 - exp(-x_vec3/2) ./ (2*pi) .* I_int);
    end
    fprintf('%-30s  %10.4f  %10.4f  %10.4f\n', configs{k,1}, C_fft3);
end
print_EX_3(x_vec3, C_quad3, C_res3, C_mc3, C_fft3);

%% Ex 4 - Normal Mean-Variance Mixture (NIG, alpha = 1/2)
fprintf('===== Exercise 4 =====\n');
sigma = 0.20; kappa = 1; eta = 3; alpha = 1/2; M = 15; N4 = 2^M;
x_vec4 = (-0.25:0.01:0.25).';
phi_ex4 = charFun_NMVM(sigma, kappa, eta, T, alpha);

% (a) FFT - Choice 1: small moneyness step (dz = 0.0025)
dz1 = 0.0025; dx1 = 2*pi / (N4*dz1);
[xg41, If41] = compute_FFT(phi_ex4, N4, dx1);
C_fft41 = B * F0 * (1 - exp(-x_vec4/2) ./ (2*pi) .* interp1(xg41, If41, x_vec4, 'spline'));

% (a) FFT - Choice 2: extremum in frequency domain (u_max = 500)
dx2 = 1000 / N4;
[xg42, If42] = compute_FFT(phi_ex4, N4, dx2);
C_fft42 = B * F0 * (1 - exp(-x_vec4/2) ./ (2*pi) .* interp1(xg42, If42, x_vec4, 'spline'));

% (a) FFT - Choice 3: perfectly aligned with the 1% moneyness grid (dz = 0.01)
dz3 = 0.01; dx3 = 2*pi / (N4*dz3);
[xg43, If43] = compute_FFT(phi_ex4, N4, dx3);
C_fft43 = B * F0 * (1 - exp(-x_vec4/2) ./ (2*pi) .* interp1(xg43, If43, x_vec4, 'spline'));

% (b) Quadrature & (c) Monte Carlo
C_quad4 = Quadrature(x_vec4, phi_ex4, B, F0);
C_mc4 = priceMC_ex4(x_vec4, F0, B, sigma, kappa, eta, T, 1e6, alpha);

% Outputs: Plots and errors
fig = comparison_plot_ex4(x_vec4, C_quad4, C_fft41, C_fft42, C_fft43, C_mc4);
[err1, err2, err3, errMC] = errors_ex4(C_quad4, C_fft41, C_fft42, C_fft43, C_mc4);