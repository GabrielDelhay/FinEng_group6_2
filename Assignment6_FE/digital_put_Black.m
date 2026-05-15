function p = digital_put_Black(F, K, sigma, tau)
% Risk-neutral digital put probability under Black (1976) lognormal dynamics
% INPUTS:
    % F     : forward price F(t0, T_reset)
    % K     : strike price
    % sigma : implied Black volatility interpolated at K
    % tau   : time to reset date in years
%
% OUTPUT:
    % p     : risk-neutral probability P(S(tau) < K) = N(-d2)
    d2 = (log(F/K) - 0.5*sigma^2*tau) / (sigma*sqrt(tau));
    p  = normcdf(-d2);
end