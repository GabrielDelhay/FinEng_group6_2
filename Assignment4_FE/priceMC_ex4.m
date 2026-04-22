function C = priceMC_ex4(x, F0, B, sigma, kappa, eta, T, N, alpha)
% Exact MC for NMVM alpha = 1/2 (NIG).
%   Z = sigma^2 * G_T ~ InverseGaussian(T*sigma^2, T^2*sigma^2/kappa)
%   f_T = -T*lnL[eta] + sqrt(Z)*W - (eta + 1/2)*Z

rng(42);

% martingale drift
lnL_eta = (T/kappa) * ((1-alpha)/alpha) * ...
          (1 - (1 + eta*kappa*sigma^2/(1-alpha))^alpha);

% IG parameters matched to Baviera's subordinator via the Laplace transform
m = T * sigma^2;
l = T^2 * sigma^2 / kappa;
Z = random(makedist('InverseGaussian','mu',m,'lambda',l), N, 1);
W = randn(N, 1);

% log-return at maturity
fT = -T*lnL_eta + sqrt(Z).*W - (eta + 0.5).*Z;

% vectorised over x (same trick as ex3)
x_row  = x(:).';                          % force row for broadcasting
payoff = min(exp(fT), exp(-x_row));       % Nx1 vs 1xM -> NxM
C      = B * F0 * (1 - mean(payoff, 1));  % 1xM
C      = reshape(C, size(x));             % preserve input shape

end