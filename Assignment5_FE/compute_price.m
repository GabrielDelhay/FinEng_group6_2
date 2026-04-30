function NPV_tot = compute_price(N, d, disc, t0, strikes, maturities, flat_vols, bond, spread)

    [spot_vols, ~, B_cap, fwd_rates, delta_fwd, tau_expiry, ~, cap_idx] = ...
        lmm_spot_vols(flat_vols, strikes, maturities, d, disc, t0);

    % Reconstruct schedule
    i_10y = cap_idx(10) + 1;

    s.B_pay     = B_cap(1 : i_10y);
    s.delta_pay = delta_fwd(1 : i_10y);
    s.fwd_rates = fwd_rates(1 : i_10y);
    s.tau_reset = tau_expiry(1 : i_10y);
    s.i_first   = 2;
    s.i_3y      = cap_idx(3)  + 1;
    s.i_6y      = cap_idx(6)  + 1;
    s.i_10y     = i_10y;

    m.spot_vols   = spot_vols(1 : i_10y, :);
    m.strikes_mkt = strikes / 100;

    % Compute final NPV
    NPV_A   = price_floating_leg(N, s.B_pay, s.delta_pay, spread);
    NPV_B   = price_coupon_leg(N, s, m, bond);
    NPV_tot = NPV_A - NPV_B;
end