function lnL = laplace_exponent(omega, p, alpha, T)
    % Laplace exponent of the tempered stable subordinator G 
    % ln L[omega] = (T/kappa) * (1-alpha)/alpha * [1 - (1 + omega*kappa*sigma^2/(1-alpha))^alpha]
    %
    % Inputs:
    %   omega : can be a scalar or vector (real or complex)
    %   p     : [sigma, kappa, eta]
    %   alpha : stability index (2/3 in our case)
    %   T     : time to maturity (= Delta_t in slide notation)
    
    sigma = p(1);
    kappa = p(2);
    
    lnL = (T / kappa) .* (1 - alpha) / alpha .* ...
          (1 - (1 + omega .* kappa .* sigma^2 ./ (1 - alpha)).^alpha);
end
