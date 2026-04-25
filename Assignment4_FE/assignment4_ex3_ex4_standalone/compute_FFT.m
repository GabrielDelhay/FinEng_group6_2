function [x, I_fft] = compute_FFT(phi, N, du)
% compute_FFT - Lewis integral I(x) on a moneyness grid via FFT.
%
% Inputs:
%   phi - characteristic function handle (martingale measure, phi(-i)=1)
%   N   - number of FFT nodes (power of two)
%   du  - frequency step. The moneyness step is dx = 2*pi/(N*du).
%
% Outputs:
%   x     - moneyness grid (Nx1)
%   I_fft - integral I(x) = integral over R of g(u) du, where g is the
%           Lewis integrand. The call price is then
%           C(x) = B*F0*(1 - exp(-x/2)/(2*pi).*I_fft(x)).

dz = 2*pi / (N*du);
u  = (0:N-1).' * du;

b  = N * dz / 2;
x  = -b + (0:N-1).' * dz;

% integrand without the Fourier kernel (see lewisIntegrand): the kernel
% exp(-i*x*u) is supplied by the FFT itself.
psi = phi(-u - 1i/2) ./ (u.^2 + 0.25);

% Simpson weights
w           = ones(N,1);
w(2:2:end)  = 4;
w(3:2:end)  = 2;
w(1)        = 1;
w(end)      = 1;
w           = w * du / 3;

fft_input  = psi .* w .* exp(-1i * u * x(1));
fft_output = real(fft(fft_input));

% the integrand is Hermitian => integral on R = 2 * Re(integral on [0,inf))
I_fft = 2 * fft_output;

end