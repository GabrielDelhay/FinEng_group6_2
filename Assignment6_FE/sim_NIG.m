function X = sim_NIG(N_sim, p, alpha, tau)
% Simulates N_sim paths of the NIG log-return X(tau) = ln(S(tau)/F0)
% under the risk-neutral measure, with martingale correction ensuring E[e^X] = 1.
% Uses the NIG representation as a Brownian motion time-changed by an Inverse Gaussian.
%
% INPUTS:
    % N_sim  : number of Monte Carlo paths
    % p      : calibrated NIG parameters [sigma, kappa, eta]
    % alpha  : mixing exponent = 0.5 (NIG), passed through to laplace_exponent
    % tau    : time horizon in years
% OUTPUT:
    % X      : (N_sim x 1) vector of simulated log-returns ln(S(tau)/F0)
%    
    sigma = p(1);  kappa = p(2);  eta = p(3);

    % G ~ IG(mu = tau*sigma^2, lambda = tau^2*sigma^2/kappa)
    mu_G  = tau * sigma^2;
    lam_G = tau^2 * sigma^2 / kappa;

    % Simulate Inverse Gaussian
    v  = randn(N_sim, 1);
    y  = v.^2;
    x  = mu_G + mu_G^2*y/(2*lam_G) - mu_G/(2*lam_G)*sqrt(4*mu_G*lam_G*y + mu_G^2*y.^2);
    u  = rand(N_sim, 1);
    G  = x.*(u <= mu_G./(mu_G+x)) + mu_G^2./x.*(u > mu_G./(mu_G+x));

    X_raw = -(1+2*eta)/2 * G + sqrt(G) .* randn(N_sim, 1);

    % Martingale correction: X = -lnL(eta) + X_raw
    lnL_eta = laplace_exponent(eta, p, alpha, tau);
    X = -lnL_eta + X_raw;
end