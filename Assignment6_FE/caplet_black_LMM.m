function price = caplet_black_LMM(L, K, delta, B_pay, tau, sigma)
% CAPLET_BLACK_LMM  Vectorized Black caplet pricer.
%
%   All inputs may be scalars or arrays of compatible size (MATLAB
%   element-wise broadcasting). The output has the broadcast shape.
%
% INPUTS:
%   L      : forward rate(s) of each caplet
%   K      : strike(s)
%   delta  : Act/360 year fraction(s) of the caplet period
%   B_pay  : discount factor(s) at the payment date of each caplet
%   tau    : Act/365 time(s) from t0 to caplet reset
%   sigma  : caplet spot vol(s)
%
% OUTPUT:
%   price  : caplet price(s) per unit notional, broadcast shape


    vol_sqrt_T = sigma .* sqrt(tau);
    d1 = ( log(L ./ K) + 0.5 .* sigma.^2 .* tau ) ./ vol_sqrt_T;
    d2 = d1 - vol_sqrt_T;

    price = B_pay .* delta .* ( L .* normcdf(d1) - K .* normcdf(d2) );
end