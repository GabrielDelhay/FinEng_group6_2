function D = stochasticDiscount_fwd(x, deltaX, a, sigma, dt, D_curr)
% Realized stochastic one-step discount per (node, branch) - forward discretization.
%
% Implements the slide formula (forward discretization):
%
%   D / D_curr = exp{ -0.5*(sigma_star)^2 - (sigma_star/sigma_hat) * (deltaX + mu_hat*x) }
%
% where:
%   D_curr     = D_{t_i}(t_{i+1}; x_i): stochastic ZCB at the current node,
%                from HJM_ZCB_vec — same for all branches.
%   D          = D(t_i, t_{i+1}; x_i, branch_k): realized stochastic discount
%                specific to each branch (depends on Delta x_{i+1}).
%
% INPUTS
%   x       - (N x 1) OU state at t_i  (column)
%   deltaX  - (N x 3) realized jump x_{i+1} - x_i on each branch
%             columns ordered as [up_target, mid_target, down_target] to match
%             p_inner / p_bot / p_top from treeProbabilities.m
%   a       - HW mean reversion
%   sigma   - HW volatility
%   dt      - tree time step (in years)
%   D_curr  - (N x 1) one-period stochastic ZCB at the node,
%             precomputed via HJM_ZCB_vec(x, t_i, t_{i+1}, B0_t, B0_tnext, a, sigma)
%
% OUTPUT
%   D       - (N x 3) realized stochastic discount per (node, branch)

x      = x(:);
D_curr = D_curr(:);

%% Step-constant coefficients
mu_hat   = 1 - exp(-a*dt);
sig_hat  = sigma * sqrt( (1 - exp(-2*a*dt)) / (2*a) );
sig_star = (sigma/a) * sqrt( dt ...
    - 2*(1 - exp(-a*dt))/a ...
    + (1 - exp(-2*a*dt))/(2*a) );

%% Broadcasting: x (N x 1), deltaX (N x 3) -> N x 3
expo = -0.5*sig_star^2 ...
    - (sig_star/sig_hat) .* ( deltaX + mu_hat .* x );

D = D_curr .* exp(expo);     % (N x 1) .* (N x 3) -> (N x 3)
end
