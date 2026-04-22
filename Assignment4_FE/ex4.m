clear; clc; close all;

[F0, B, r, T] = get_datas();

% Import the info from the assignment
sigma = 0.20;  kappa = 1;  eta = 3;  alpha = 1/3; % Note: Ex 4d uses 1/3, Ex 4 main uses 1/2
x_vec = (-0.25:0.01:0.25)';

% FFT Base Parameters
M = 15; 
N = 2^M;

%% Choice 1: Small moneyness step (dz = 0.0025)
dz1 = 0.0025; 
dx1 = 2*pi / (N * dz1); 
[x_grid1, I_fft1] = compute_FFT_ex4(sigma, kappa, eta, T, alpha, N, dx1);
I_interp1 = interp1(x_grid1, I_fft1, x_vec, 'spline');
C_fft1 = B * F0 * (1 - exp(-x_vec/2) ./ (2*pi) .* I_interp1);

%% Choice 2: Extremum in frequency domain (x_max = 500)
% Total domain length = 1000 -> N * dx = 1000
dx2 = 1000 / N;  
dz2 = 2*pi / (N * dx2);
[x_grid2, I_fft2] = compute_FFT_ex4(sigma, kappa, eta, T, alpha, N, dx2);
I_interp2 = interp1(x_grid2, I_fft2, x_vec, 'spline');
C_fft2 = B * F0 * (1 - exp(-x_vec/2) ./ (2*pi) .* I_interp2);

%% Choice 3: Optimal (Exact 1% moneyness grid, dz = 0.01)
dz3 = 0.01;
dx3 = 2*pi / (N * dz3);
[x_grid3, I_fft3] = compute_FFT_ex4(sigma, kappa, eta, T, alpha, N, dx3);
% Direct subset matching or interpolation (safe fallback)
I_interp3 = interp1(x_grid3, I_fft3, x_vec, 'spline');
C_fft3 = B * F0 * (1 - exp(-x_vec/2) ./ (2*pi) .* I_interp3);

%% Benchmark calculations
% Quadrature method
C_quad = Quadrature_ex4(x_vec, sigma, kappa, eta, T, alpha, B, F0);
% Monte Carlo method
C_mc = priceMC_ex4(x_vec, F0, B, sigma, kappa, eta, T, 1e6, alpha);

%% Plots & Comparisons
figure('Position', [100, 100, 800, 500]);
plot(x_vec*100, C_quad, 'k-', 'LineWidth', 2, 'DisplayName','Quadrature (Exact)'); hold on;
plot(x_vec*100, C_fft1, 'bo', 'MarkerSize', 6, 'DisplayName','FFT Choice 1 (dz=0.0025)');
plot(x_vec*100, C_fft2, 'r+', 'MarkerSize', 6, 'DisplayName','FFT Choice 2 (Freq Extr=500)');
plot(x_vec*100, C_fft3, 'g.', 'MarkerSize', 12, 'DisplayName','FFT Choice 3 (dz=0.01)');

xlabel('Moneyness x = log(F_0/K) [%]'); 
ylabel('Call Price');
title('Call Prices Comparison: FFT Choices vs Quadrature');
legend('Location','northwest'); 
grid on;

% Error analysis
err1 = max(abs(C_quad - C_fft1));
err2 = max(abs(C_quad - C_fft2));
err3 = max(abs(C_quad - C_fft3));
fprintf('Max Error Choice 1: %e\n', err1);
fprintf('Max Error Choice 2: %e\n', err2);
fprintf('Max Error Choice 3: %e\n', err3);