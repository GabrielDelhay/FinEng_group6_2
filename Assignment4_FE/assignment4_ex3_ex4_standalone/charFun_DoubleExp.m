function phi = charFun_DoubleExp(p_plus, p_minus)
% charFun_DoubleExp - Characteristic function used in Exercise 3.
%
%   phi^c(u) = exp(i*mu*u) / ((1 - i*u/p_plus) * (1 + i*u/p_minus))
%
% mu is fixed by the martingale condition phi^c(-i) = 1, which yields
%   mu = log( (1 - 1/p_plus) * (1 + 1/p_minus) ).
%
% Inputs:
%   p_plus  - right tail parameter (must satisfy p_plus > 1)
%   p_minus - left  tail parameter
%
% Output:
%   phi - function handle phi(v), with v real or complex.

mu  = log((1 - 1/p_plus) * (1 + 1/p_minus));
phi = @(v) exp(1i * mu * v) ./ ...
          ((1 - 1i .* v ./ p_plus) .* (1 + 1i .* v ./ p_minus));

end
