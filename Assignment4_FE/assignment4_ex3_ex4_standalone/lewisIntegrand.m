function g = lewisIntegrand(u, x, phi)
% lewisIntegrand - Generic integrand of the Lewis formula.
%
% The integrand g(u) = phi(-u - i/2) * exp(-i*x*u) / (u^2 + 1/4)
% is shared by every Lewis-type pricing technique (quadrature, FFT, ...).
% The model dependence is fully encapsulated in the characteristic
% function handle phi, so that swapping model only requires building a
% different phi outside this routine.
%
% Inputs:
%   u   - integration variable (real vector)
%   x   - moneyness log(F0/K)
%   phi - function handle: characteristic function of f_T = ln(F_T/F_0)
%         under the martingale measure. Must satisfy phi(-i) = 1.
%         Called as phi(v) with v complex.
%
% Output:
%   g   - integrand value (complex). Take real part outside if needed.

v = -u - 1i/2;
g = phi(v) .* exp(-1i .* x .* u) ./ (u.^2 + 0.25);

end
