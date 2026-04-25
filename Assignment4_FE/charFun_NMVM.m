function phi = charFun_NMVM(sigma, kappa, eta, T, alpha)
% charFun_NMVM - Characteristic function of the Normal Mean-Variance
% Mixture used in Exercise 4.
%
%   phi(v) = exp( -i*v*lnL[eta] + lnL[ (v^2 + i*(1+2*eta)*v)/2 ] )
%
% with the Laplace exponent of a tempered-stable subordinator
%   lnL[w] = (T/kappa) * (1-alpha)/alpha * [1 - (1 + w*kappa*sigma^2/(1-alpha))^alpha].
%
% The first term enforces the martingale condition phi(-i) = 1.
%
% Inputs:
%   sigma, kappa, eta - NMVM parameters
%   T                 - time to maturity
%   alpha             - stability index (1/2 for NIG, 1/3 for the variant)
%
% Output:
%   phi - function handle phi(v).

lnL = @(w) (T/kappa) .* ((1-alpha)/alpha) .* ...
          (1 - (1 + w .* kappa * sigma^2 / (1-alpha)).^alpha);

lnL_eta = lnL(eta);

phi = @(v) exp(-1i .* v .* lnL_eta + ...
               lnL((v.^2 + 1i .* (1 + 2*eta) .* v) / 2));

end
