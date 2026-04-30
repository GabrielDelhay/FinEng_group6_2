function price = digital_black(L, K, delta, B_pay, tau, sigma)
%  INPUTS:
%    L      : forward rate L_i(t0)
%    K      : strike
%    delta  : Act/360 year fraction of the coupon period
%    B_pay  : discount factor to payment date B(0, T_{i+1})
%    tau    : Act/365 time to expiry  (T_i - t0)/365
%    sigma  : spot vol of this caplet (Black, on L)
%
%  OUTPUT:
%    price  : digital price per unit notional, in "rate units"

vol_sqrt_T = sigma .* sqrt(tau);
d2 = (log(L./K) - 0.5 .* sigma.^2 .* tau) ./ vol_sqrt_T;

price = B_pay .* delta .* normcdf(d2);
end