function [C, SE, CI] = priceBasketCall_MC(S1_0, S2_0, sigma1, sigma2, d1, d2, rho, r, T, M)
% priceBasketCall_MC  Monte Carlo price of a European call on an equally weighted basket of two lognormal stocks (ENI, AXA)
%
% INPUTS:
%   S1_0, S2_0     : spot prices at t0       (ENI, AXA)
%   sigma1, sigma2 : Black volatilities      (annualized)
%   d1, d2         : continuous dividend yields
%   rho            : correlation between the two log-returns
%   r              : continuous-compounding zero rate, r = -log(B(t0,T))/T
%   T              : time to maturity in years
%   M              : number of Monte Carlo paths
%
% OUTPUTS:
%   C   : MC estimate of the basket call price (per unit of notional)
%   SE  : MC standard error
%   CI  : 95% confidence interval [lower, upper]

%% Discount factor consistent with the zero rate and two correlated standard normals  (1 x M each)
B = exp(-r * T);

Z1 = randn(1, M);
W  = randn(1, M);
Z2 = rho * Z1 + sqrt(1 - rho^2) * W;

%% Risk-neutral simulation of the ratios S_n(T) / S_n(0)   (slide 6)
ratio1 = exp((r - d1 - 0.5*sigma1^2)*T + sigma1*sqrt(T)*Z1);
ratio2 = exp((r - d2 - 0.5*sigma2^2)*T + sigma2*sqrt(T)*Z2);

%% Equally weighted normalized basket and discounted payoff
Bsk = 0.5 * (ratio1 + ratio2);        % 1 x M
payoff = max(0, Bsk - 0.95);          % 1 x M
discPay = B * payoff;                 % 1 x M

%% Monte Carlo statistics
C  = mean(discPay);
SE = std(discPay) / sqrt(M);
CI = C + 1.96 * SE * [-1, 1];

end 