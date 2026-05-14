function [D_t] = HJM_ZCB_vec(x, t, T, B0_t, B0_T, a, sigma)
% HJM_ZCB_vec  Vectorized stochastic ZCB prices in 1-factor Hull-White.
%
% Computes D_t(T;x) = the stochastic ZCB observed at time t, for an
% array of maturities T, evaluated at every state-variable node x.
% Uses the closed-form HJM Gaussian solution:
%
%     D_t(T;x) = B_t(T) * exp( - sigma_hat(t,T)*x - 0.5 * V(t,T) )
%
% where:
%   B_t(T) = B0(T)/B0(t)                              (deterministic forward DF)
%   sigma_hat(t,T) = (1 - exp(-a*(T-t)))/a            (HW state coefficient)
%   V(t,T) = closed-form HJM variance integral        (convexity correction)
%
% NOTATION
%   B0(.)    : initial discount factor curve (ZCB price observed at t=0)
%   B_t(T)   : forward discount factor for [t,T], deterministic from B0
%   D_t(T;x) : stochastic ZCB at tree node x, observed at t, maturity T
%
% PARAMETERS
%   x       - state variables at time t (N elements, forced to column)
%   t       - observation time in years (scalar)
%   T       - vector of maturities in years (M elements)
%   B0_t    - initial discount factor B0(t) (scalar)
%   B0_T    - initial discount factors B0(T_k) (M elements)
%   a       - HW mean reversion rate
%   sigma   - HW volatility
%
% RETURNS
%   D_t     - matrix of stochastic ZCB prices, size N x M

    %% 1. Array reshaping
    x    = x(:);                            % column N x 1
    T    = T(:)';                           % row    1 x M
    B0_T = B0_T(:)';                        % row    1 x M

    %% 2. Analytical pieces
    % Closed-form variance integral V(t,T)
    term1 = (2*(1 - exp(-a*t))/a)    .* (1 - exp(-a*(T - t)));
    term2 = ((1 - exp(-2*a*t))/(2*a)).* (1 - exp(-2*a*(T - t)));
    V_tT  = (sigma^2 / a^2) * (term1 - term2);

    % State-loading coefficient sigma_hat(t,T)
    sigma_hat = (1 - exp(-a*(T - t))) / a;

    % Forward discount factor B_t(T) = B0(T)/B0(t) (deterministic)
    B_t = B0_T ./ B0_t;

    %% 3. Stochastic ZCB matrix (N x M via broadcasting)
    D_t = B_t .* exp(-x .* sigma_hat - 0.5 * V_tT);

    %% 4. Boundary: D_t(t;x) = 1 for any x when T equals t
    D_t(:, abs(T - t) < 1e-12) = 1;
end
