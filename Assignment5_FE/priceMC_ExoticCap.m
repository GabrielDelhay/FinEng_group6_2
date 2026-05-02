function payoffs = priceMC_ExoticCap(nu, B0_T, delta, L0, C, RHO, N, Nsim)
% Simulates the exotic cap price under the spot measure via Monte Carlo.
% Forward bonds follow log-normal dynamics with BMM spot vols nu and
% exponential correlation structure encoded in C (Cholesky of RHO).

% Inputs: nu (BMM spot vols, 1x15)
%         B0_T (discount factors, 1x17),
%         delta (year fractions, 1x16), 
%         L0 (forward Libors, 1x16),
%         C (Cholesky lower factor, 16x16), 
%         RHO (correlation matrix, 16x16),
%         N (number of periods), 
%         Nsim (number of Monte Carlo paths) 

    nu_col      = nu(:);
    B_sim       = repmat((B0_T(2:end) ./ B0_T(1:end-1))', 1, Nsim);
    fixedLibors = zeros(N, Nsim);
    fixedLibors(1, :) = L0(1);
    cumDiscount = ones(1, Nsim);
    payoffs     = zeros(1, Nsim);

    for k = 0:N-1
        dt    = delta(k+1);
        dW    = sqrt(dt) * (C * randn(N, Nsim));
        B_old = B_sim;
        nr    = (k+1):(N-1);

        if ~isempty(nr)
            nu_a           = nu_col(nr);
            drift_a        = -(nu_a .* (tril(RHO(nr,nr), -1) * nu_a));
            B_sim(nr+1, :) = B_old(nr+1, :) .* exp((drift_a - 0.5*nu_a.^2)*dt - nu_a .* dW(nr, :));
        end

        if k <= N-2
            fixedLibors(k+2, :) = (1./B_sim(k+2, :) - 1) / delta(k+2);
        end

        cumDiscount = cumDiscount .* B_old(k+1, :);

        if k >= 1 && k <= N-2
            payoff  = delta(k+1) * max(fixedLibors(k+1,:) - fixedLibors(k,:) - 0.0005, 0.0);
            payoffs = payoffs + cumDiscount .* B_sim(k+2, :) .* payoff;
        end
    end
    payoffs = payoffs';
end
