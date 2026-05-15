function lnL = laplace_exponent(omega, p, alpha, T)
% Laplace exponent of the NMV subordinator G, used for martingale correction
% and characteristic function computation.
%
% INPUTS:
    % omega  : Laplace argument (scalar or vector, real or complex)
    % p      : NIG parameters [sigma, kappa, eta]
    % alpha  : mixing exponent (= 0.5 for NIG, = 2/3 in general NMV)
    % T      : time horizon in years
%
% OUTPUT:
    % lnL    : ln E[e^{-omega*G}] = (T/kappa) * (1-alpha)/alpha
    %          * [1 - (1 + omega*kappa*sigma^2/(1-alpha))^alpha]
    
    sigma = p(1);
    kappa = p(2);
    
    lnL = (T / kappa) .* (1 - alpha) / alpha .* ...
          (1 - (1 + omega .* kappa .* sigma^2 ./ (1 - alpha)).^alpha);
end
