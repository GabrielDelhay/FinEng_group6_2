function [NPV_A, NPV_B_flat, X_flat, NPV_B_smile, X_smile] = price_structured_bond(N, spread, bond, B_cap, delta_fwd, tau_expiry, fwd_rates, cap_maturity_idx, spot_vols, strikes)
% PRICE_STRUCTURED_BOND Computes the NPVs and upfronts for the structured bond

    % 1. Build payment schedule and regime boundaries
    i_3y  = cap_maturity_idx(3)  + 1;     % Period paying at 3y
    i_6y  = cap_maturity_idx(6)  + 1;     % Period paying at 6y
    i_10y = cap_maturity_idx(10) + 1;     % Period paying at 10y

    schedule.B_pay     = B_cap(1 : i_10y);
    schedule.delta_pay = delta_fwd(1 : i_10y);
    schedule.tau_reset = tau_expiry(1 : i_10y);
    schedule.fwd_rates = fwd_rates(1 : i_10y);
    schedule.i_first   = 2;               % First structured coupon
    schedule.i_3y      = i_3y;
    schedule.i_6y      = i_6y;
    schedule.i_10y     = i_10y;

    % 2. Setup market data
    market.spot_vols   = spot_vols(1 : i_10y, :);
    market.strikes_mkt = strikes / 100;

    % 3. Price floating leg (Party A)
    NPV_A = price_floating_leg(N, schedule.B_pay, schedule.delta_pay, spread);

    % 4. Price coupon leg (Party B) - Flat Black
    bond.use_smile = false;
    NPV_B_flat = price_coupon_leg(N, schedule, market, bond);
    X_flat     = (NPV_A - NPV_B_flat) / N;

    % 5. Price coupon leg (Party B) - Smile Corrected
    bond.use_smile = true;
    NPV_B_smile = price_coupon_leg(N, schedule, market, bond);
    X_smile     = (NPV_A - NPV_B_smile) / N;

end