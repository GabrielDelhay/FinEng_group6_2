function g = lewisIntegrand(u, x, p_plus, p_minus, mu)
%% lewisIntegrand  -  Integrand of the Lewis formula
%  Inputs:
%    u      - integration variable (real vector)
%    x      - moneyness = log(F/K)
%    p_plus - right tail parameter  (p+ > 1 required)
%    p_minus- left tail parameter
%    mu     - martingale drift
%
%  Output:
%    val    - real part of integrand (scalar or vector, same size as u)

v = - u - 1i/2;
phi = exp(1i * mu * v) ./ ((1 - 1i * v / p_plus) .* (1 + 1i * v / p_minus));

g = phi .* exp(- 1i * x * u) ./ (u.^2 + 0.25);

end
