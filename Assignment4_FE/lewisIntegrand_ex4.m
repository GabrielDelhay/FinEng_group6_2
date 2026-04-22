function g = lewisIntegrand_ex4(u, x, sigma, kappa, eta, T, alpha)
% Lewis integrand for the NMVM model (ex4).
%   phi(v) = exp(-i*v*lnL[eta]) * L[(v^2 + i(1+2eta)v)/2]
%   ln L[w] = (T/k)*(1-a)/a * (1 - (1 + w*k*sigma^2/(1-a))^a)
% Martingale drift is embedded via lnL[eta] (no external mu needed).

v = -u - 1i/2;

lnL_eta = (T/kappa) * ((1-alpha)/alpha) * ...
          (1 - (1 + eta*kappa*sigma^2/(1-alpha))^alpha);

w     = (v.^2 + 1i*(1 + 2*eta).*v) / 2;
lnL_w = (T/kappa) * ((1-alpha)/alpha) * ...
        (1 - (1 + w*kappa*sigma^2/(1-alpha)).^alpha);

phi = exp(-1i * v * lnL_eta + lnL_w);
g   = phi .* exp(-1i * x * u) ./ (u.^2 + 0.25);

end
