function X = sim_NIG(N_sim, p, alpha, tau)
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

    % X_raw = -(1+2*eta)/2 * G + sqrt(G) * W  (coeff BM = 1, pas sigma)
    X_raw = -(1+2*eta)/2 * G + sqrt(G) .* randn(N_sim, 1);

    % Martingale correction: X = -lnL(eta) + X_raw
    lnL_eta = laplace_exponent(eta, p, alpha, tau);   % < 0 pour eta > 0
    X = -lnL_eta + X_raw;
end