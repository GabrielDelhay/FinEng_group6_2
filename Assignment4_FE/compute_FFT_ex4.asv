function [x, C_fft] = compute_FFT_ex4(sigma, kappa, eta, T, alpha, N, du)

lambda = 2*pi / (N*du);

% frequency grid
u = (0:N-1)' * du;

% moneyness grid
b = N * lambda / 2;
x = -b + (0:N-1)' * lambda;

% phi(v) at v = -u - i/2   (martingale drift embedded via lnL[eta])
v = -u - 1i/2;

lnL_eta = (T/kappa) * ((1-alpha)/alpha) * ...
          (1 - (1 + eta*kappa*sigma^2/(1-alpha))^alpha);

z     = (v.^2 + 1i*(1 + 2*eta).*v) / 2;
lnL_z = (T/kappa) * ((1-alpha)/alpha) * ...
        (1 - (1 + z*kappa*sigma^2/(1-alpha)).^alpha);

phi = exp(-1i * v * lnL_eta + lnL_z);

% FFT integrand
psi = phi ./ (u.^2 + 0.25);

% Simpson weights
w = ones(N,1);
w(2:2:end) = 4;
w(3:2:end) = 2;
w(1)   = 1;
w(end) = 1;
w = w * du / 3;

% FFT
fft_input  = psi .* w .* exp(-1i * u * x(1));
fft_output = real(fft(fft_input));

% integrand is Hermitian -> integral on R = 2 * Re(integral on [0,inf))
C_fft = 2 * fft_output;

end
