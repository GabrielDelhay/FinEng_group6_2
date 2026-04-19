function C = priceMC(x, F0, B, p_plus, p_minus, mu, N)
%% priceMC  -  Monte Carlo pricing via exact simulation
%
%  Inputs:
%    x_vec   - moneyness values log(F/K)  [vector]
%    F, B    - forward, discount factor
%    p_plus, p_minus, mu - model parameters
%    N_paths - number of MC paths (e.g. 1e6)
%
%  Outputs:
%    C    - call prices  [same size as x_vec]

rng(42);  % reproducibility

% Simulation de fT = mu + E+ - E-
E_plus  = exprnd(1/p_plus,  N, 1);
E_minus = exprnd(1/p_minus, N, 1);
fT       = mu + E_plus - E_minus;

payoff = min(exp(fT), exp(-x));

% Prix = B * E[payoff]
C = B * F0 * (1 - mean(payoff));

end
