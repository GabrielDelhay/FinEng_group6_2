function [x, C_fft] = compute_FFT(p_plus, p_minus, mu)

i = 1i;

% paramètres FFT
N   = 2^12;        % nombre de points
eta = 0.05;        % pas en u
lambda = 2*pi/(N*eta);

% grille en u
u = (0:N-1)' * eta;

% grille en x
b = N * lambda / 2;
x = -b + (0:N-1)' * lambda;

% fonction phi(v)
v = -u - 1i/2;

phi = exp(1i * mu * v) ./ ...
    ((1 - 1i * v / p_plus) .* (1 + 1i * v / p_minus));

% intégrande FFT
psi = phi ./ (u.^2 + 0.25);

% poids de Simpson (important !)
w = ones(N,1);
w(2:2:end) = 4;
w(3:2:end) = 2;
w(1) = 1;
w(end) = 1;
w = w * eta / 3;

% FFT
fft_input = psi .* w .* exp(-1i * u * x(1));
fft_output = real(fft(fft_input));

C_fft = fft_output;

end