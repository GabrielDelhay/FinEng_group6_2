function NPV_B = price_coupon_leg(notional, schedule, market, bond)
%  COUPON STRUCTURE:
%    i = 1            : 4% fixed
%    i = 2..i_3y      : Regime 1: L + 1.00% if L<=4.20% else 4.50%
%                       = (L+1%) - (L-4.20%)+ - 0.70% * 1{L>4.20%}
%    i_3y+1..i_6y     : Regime 2: L + 1.20% if L<=4.70% else 4.90%
%                       = (L+1.20%) - (L-4.70%)+ - 0.70% * 1{L>4.70%}
%    i_6y+1..i_10y    : Regime 3, two variants in the termsheet:
%        (3a) L + 1.10%  CAPPED at 5.10%
%             = (L+1.10%) - (L-4.00%)+
%        (3b) L + 1.30%  if L<=5.40% else 5.60%
%             = (L+1.30%) - (L-5.40%)+ - 0.90% * 1{L>5.40%}
%
%  INPUTS:
%    notional  : N 
%
%    schedule  : struct with fields (all column vectors of length N_total
%                aligned with EX_1.m indexing on periods)
%       .B_pay      : (N x 1) discount factors B(t0, T_i)  for i=1..N
%       .delta_pay  : (N x 1) Act/360 year fractions delta_i
%       .tau_reset  : (N x 1) Act/365 times to reset (T_{i-1}-t0)/365
%       .fwd_rates  : (N x 1) forward rates L_{i-1}(t0)
%       .i_first    : index of the first STRUCTURED coupon
%       .i_3y       : index of the last period in regime 1 (<= 3y)
%       .i_6y       : index of the last period in regime 2 (<= 6y)
%       .i_10y      : index of the last period (= 10y, end of bond)
%
%    market    : struct
%       .spot_vols    : (N x N_strikes) calibrated LMM spot-vol surface
%       .strikes_mkt  : (1 x N_strikes) market strikes [decimal]
%
%    bond      : struct
%       .first_cpn      : value of the first fixed coupon 
%       .spread1, .K1, .c_dig1 : Regime 1 parameters
%       .spread2, .K2, .c_dig2 : Regime 2 parameters
%       .regime3_type   : '3a' or '3b'
%       .spread3, .K3, .c_dig3 : Regime 3 parameters (c_dig3 = 0 if 3a)
%
%  OUTPUT:
%    NPV        : present value of the coupon leg, in currency units

%% Pre-interpolate spot vols at the relevant strikes.
% We interpolate ONCE per strike (vectorised across all caplets) and then index the resulting column by the period i

sigma_K1 = interp1(market.strikes_mkt(:), market.spot_vols.', bond.K1, 'spline').';
sigma_K2 = interp1(market.strikes_mkt(:), market.spot_vols.', bond.K2, 'spline').';
sigma_K3 = interp1(market.strikes_mkt(:), market.spot_vols.', bond.K3, 'spline').';

%% First fixed coupon (4%)

pv_first = notional * schedule.delta_pay(1) * bond.first_cpn * schedule.B_pay(1);

%% Regime 1. i = i_first : i_3y

pv_reg1 = 0;
for i = schedule.i_first : schedule.i_3y
    L = schedule.fwd_rates(i);
    delta = schedule.delta_pay(i);
    B_pay = schedule.B_pay(i);
    tau = schedule.tau_reset(i);
    sigma = sigma_K1(i);

    floater = delta * B_pay * (L + bond.spread1);
    caplet = caplet_black_LMM(L, bond.K1, delta, B_pay, tau, sigma);
    digital = bond.c_dig1 .* digital_black(L, bond.K1, delta, B_pay, tau, sigma);

    pv_reg1 = pv_reg1 + notional * (floater - caplet - digital);
end
%% Regime 2. i = i_3y+1 : i_6y

pv_reg2 = 0;
for i = schedule.i_3y + 1 : schedule.i_6y
    L = schedule.fwd_rates(i);
    delta = schedule.delta_pay(i);
    B_pay = schedule.B_pay(i);
    tau = schedule.tau_reset(i);
    sigma = sigma_K2(i);

    floater = delta * B_pay * (L + bond.spread2);
    caplet = caplet_black_LMM(L, bond.K2, delta, B_pay, tau, sigma);
    digital = bond.c_dig2 .* digital_black(L, bond.K2, delta, B_pay, tau, sigma);

    pv_reg2 = pv_reg2 + notional * (floater - caplet - digital);
end
%% Regime 3. i = i_6y+1 : i_10y

pv_reg3 = 0;
for i = schedule.i_6y + 1 : schedule.i_10y
    L = schedule.fwd_rates(i);
    delta = schedule.delta_pay(i);
    B_pay = schedule.B_pay(i);
    tau = schedule.tau_reset(i);
    sigma = sigma_K3(i);

    floater = delta * B_pay * (L + bond.spread3);
    caplet = caplet_black_LMM(L, bond.K3, delta, B_pay, tau, sigma);
    digital = bond.c_dig3 .* digital_black(L, bond.K3, delta, B_pay, tau, sigma);

    pv_reg3 = pv_reg3 + notional * (floater - caplet - digital);
end
%% NPV Party B
NPV_B = pv_first + pv_reg1 + pv_reg2 + pv_reg3;

end