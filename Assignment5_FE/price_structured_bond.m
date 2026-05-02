function [NPV_A, NPV_B_flat, X_flat, NPV_B_smile, X_smile] = price_structured_bond(N, spread, bond, B_cap, delta_fwd, tau_expiry, fwd_rates, cap_maturity_idx, spot_vols, strikes, cap_dates, t0)
% PRICE_STRUCTURED_BOND Computes the NPVs and upfronts for the structured bond

    % 1. Build payment schedule and regime boundaries
    i_3y  = cap_maturity_idx(3)  + 1;     % Period paying at 3y
    i_6y  = cap_maturity_idx(6)  + 1;     % Period paying at 6y
    i_10y = cap_maturity_idx(10) + 1;     % Period paying at 10y

    n_strikes = numel(strikes);

    schedule.B_pay     = B_cap(1 : i_10y);
    schedule.delta_pay = [yearfrac(t0, cap_dates(1), 2); delta_fwd(1 : i_10y -1)];
    schedule.tau_reset = [NaN; tau_expiry(1 : i_10y - 1)];
    schedule.fwd_rates = [NaN; fwd_rates(1 : i_10y -1)];
    schedule.i_first   = 2;               % First structured coupon
    schedule.i_3y      = i_3y;
    schedule.i_6y      = i_6y;
    schedule.i_10y     = i_10y;

    % 2. Setup market data
    market.spot_vols   = [zeros(1, n_strikes); spot_vols(1 : i_10y -1, :)];
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