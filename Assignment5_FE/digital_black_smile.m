function price = digital_black_smile(L, K, delta, B_pay, tau, sigma, dsigma_dK)
% Digital price with implied-vol smile correction.
% Replaces the flat-Black digital with the call-spread replica that accounts for the slope of the implied-vol smile at strike K:
%
%  INPUTS:
%    L         : forward rate L_i(t0)
%    K         : strike
%    delta     : Act/360 year fraction of the coupon period
%    B_pay     : discount factor to payment date B(0, T_{i+1})
%    tau       : Act/365 time to expiry  (T_i - t0)/365
%    sigma     : spot vol of this caplet at strike K (Black, on L)
%    dsigma_dK : slope of the implied-vol smile in K, evaluated at K
%
%  OUTPUT:
%    price     : digital price per unit notional, in "rate units"

vol_sqrt_T = sigma * sqrt(tau);
d1 = ( log(L/K) + 0.5 * sigma^2 * tau ) / vol_sqrt_T;
d2 = d1 - vol_sqrt_T;

price_black = B_pay * delta * normcdf(d2);
vega = B_pay * delta * L * sqrt(tau) * normpdf(d1);  % Black vega of the underlying caplet   

price = price_black - dsigma_dK * vega;  % Smile correction (call-spread slope effect)
end