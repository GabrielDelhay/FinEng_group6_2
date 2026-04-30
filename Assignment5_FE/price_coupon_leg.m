function NPV_B = price_coupon_leg(notional, schedule, market, bond)
%  COUPON STRUCTURE:
%    i = 1            : 4% fixed
%    i = 2..i_3y      : Regime 1: L + 1.00% if L<=4.20% else 4.50%
%                       = (L+1%) - (L-4.20%)+ - 0.70% * 1{L>4.20%}
%    i_3y+1..i_6y     : Regime 2: L + 1.20% if L<=4.70% else 4.90%
%                       = (L+1.20%) - (L-4.70%)+ - 0.70% * 1{L>4.70%}
%    i_6y+1..i_10y    : Regime 3: L + 1.30% if L<=5.40% else 5.60%
%                       = (L+1.30%) - (L-5.40%)+ - 0.90% * 1{L>5.40%}
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
%       .spread3, .K3, .c_dig3 : Regime 3 parameters
%       .use_smile  (opt)      : false (default) -> flat-Black digital
%                                true            -> smile-corrected digital
%
%  OUTPUT:
%    NPV_B  : present value of the coupon leg, in currency units

%% Pricing mode flag (flat Black vs smile-corrected)
if ~isfield(bond, 'use_smile') || ~bond.use_smile
    use_smile = false;
else
    use_smile = true;
end
%% Pre-interpolate spot vols at the relevant strikes.
%  Single batched call returns sigma_K of size (N x 3): column j has the caplet spot vol at strike K_j over all reset dates.
%  When use_smile is on, we additionally compute the smile slope dsigma/dK at each K_j via central differences on the same spline.

K_vec    = [bond.K1, bond.K2, bond.K3];
sigma_K  = interp1(market.strikes_mkt(:), market.spot_vols.', K_vec, 'spline').';

if use_smile
    h        = 1e-4;
    sigma_up = interp1(market.strikes_mkt(:), market.spot_vols.', K_vec+h, 'spline').';
    sigma_dn = interp1(market.strikes_mkt(:), market.spot_vols.', K_vec-h, 'spline').';
    dsigma_dk  = (sigma_up - sigma_dn) / (2*h);
else
    dsigma_dk  = zeros(size(sigma_K));   % no correction
end

%% First fixed coupon (4%)
pv_first = notional * schedule.delta_pay(1) * bond.first_cpn * schedule.B_pay(1);

%% Regime 1. i = i_first : i_3y
pv_reg1 = priceRegime(notional, schedule, schedule.i_first, schedule.i_3y, ...
                      bond.K1, bond.spread1, bond.c_dig1, ...
                      sigma_K(:,1), dsigma_dk(:,1), use_smile);

%% Regime 2. i = i_3y+1 : i_6y
pv_reg2 = priceRegime(notional, schedule, schedule.i_3y+1, schedule.i_6y, ...
                      bond.K2, bond.spread2, bond.c_dig2, ...
                      sigma_K(:,2), dsigma_dk(:,2), use_smile);

%% Regime 3. i = i_6y+1 : i_10y
pv_reg3 = priceRegime(notional, schedule, schedule.i_6y+1, schedule.i_10y, ...
                      bond.K3, bond.spread3, bond.c_dig3, ...
                      sigma_K(:,3), dsigma_dk(:,3), use_smile);

%% NPV Party B
NPV_B = pv_first + pv_reg1 + pv_reg2 + pv_reg3;

end





%% Local helper: prices one regime block

function pv = priceRegime(notional, schedule, i_start, i_end, ...
                          K_reg, spread_reg, c_dig_reg, ...
                          sigma_reg, dsigma_dk_reg, use_smile)

pv = 0;
for i = i_start : i_end
    L     = schedule.fwd_rates(i);
    delta = schedule.delta_pay(i);
    B_pay = schedule.B_pay(i);
    tau   = schedule.tau_reset(i);
    sigma = sigma_reg(i);

    floater = delta * B_pay * (L + spread_reg);
    caplet  = caplet_black_LMM(L, K_reg, delta, B_pay, tau, sigma);

    if use_smile
        digital = c_dig_reg .* digital_black_smile( ...
                       L, K_reg, delta, B_pay, tau, sigma, dsigma_dk_reg(i));
    else
        digital = c_dig_reg .* digital_black( ...
                       L, K_reg, delta, B_pay, tau, sigma);
    end

    pv = pv + notional * (floater - caplet - digital);
end
end