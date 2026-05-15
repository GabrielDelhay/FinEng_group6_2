function C = BS_call(F0, K, B, sigma, T)
% INPUTS:
    % F0    : forward price F(t0, T)
    % K     : strike price (scalar or vector)
    % B     : discount factor B(t0, T) = exp(-r*T)
    % sigma : implied Black volatility (scalar or vector, matching K)
    % T     : time to maturity in years
%
% OUTPUT:
    % C     : Black (1976) call price
    d1 = (log(F0./K) + 0.5*sigma.^2.*T) ./ (sigma.*sqrt(T));
    d2 = d1 - sigma.*sqrt(T);
    C  = B .* (F0 .* normcdf(d1) - K .* normcdf(d2));
end