function price = caplet_black_LMM(L, K, delta, B_pay, tau, sigma)
% CAPLET_BLACK_LMM  Vectorized Black caplet price under LMM.
%
%   caplet = B(0,T_{i+1}) * delta_i * [L_i * N(d1) - K * N(d2)]
%
%   with (Black 76, T-forward measure):
%     d1 = ( log(L/K) + 0.5*sigma^2*tau ) / ( sigma*sqrt(tau) )
%     d2 = d1 - sigma*sqrt(tau)
%
%   All inputs may be scalars or arrays of compatible size (MATLAB
%   element-wise broadcasting). The output has the broadcast shape.
%
% INPUTS:
%   L      : forward rate(s)         L_i(t0)
%   K      : strike(s)
%   delta  : Act/360 year fraction(s) of the caplet period
%   B_pay  : discount factor(s)      B(0, T_{i+1})
%   tau    : Act/365 time(s) to expiry, (T_i - t0)/365
%   sigma  : spot vol(s)             sigma_i
%
% OUTPUT:
%   price  : caplet price(s) per unit notional, broadcast shape

    % Numerical safeguards: avoid div-by-zero in vol_sqrt_T and log(0).
    % As sigma -> 0 the formula converges to the intrinsic value, which is
    % the correct limit of a Black-76 caplet.


    vol_sqrt_T = sigma .* sqrt(tau);
    d1 = ( log(L ./ K) + 0.5 .* sigma.^2 .* tau ) ./ vol_sqrt_T;
    d2 = d1 - vol_sqrt_T;

    price = B_pay .* delta .* ( L .* normcdf(d1) - K .* normcdf(d2) );
end