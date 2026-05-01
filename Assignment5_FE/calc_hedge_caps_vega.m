function vega_caps = calc_hedge_caps_vega(flat_vols, strikes, maturities, dates, discounts, t0, dVol, N)
% CALC_HEDGE_CAPS_VEGA Computes the 2x2 Vega matrix for 6y and 10y ATM Caps
% explicitly reusing the LMM bootstrap logic and schedule.

    % 1. Perform base bootstrap to extract the clean schedule and forward rates
    [~, ~, B_cap, fwd_rates, delta_fwd, tau_expiry, ~, cap_maturity_idx] = ...
        lmm_spot_vols(flat_vols, strikes, maturities, dates, discounts, t0);

    % Find indices in "maturities" for the 6y and 10y tenors
    idx_mat_6y = find(maturities == 6);
    idx_mat_10y = find(maturities == 10);
    
    % Find the index of the last caplet for the 6y and 10y Caps
    idx_c_6 = cap_maturity_idx(idx_mat_6y);
    idx_c_10 = cap_maturity_idx(idx_mat_10y);

    % 2. Compute ATM Strikes strictly on the LMM schedule
    % The par rate is calculated as: (1 - final_DF) / BPV
    BPV_6y = sum(delta_fwd(1:idx_c_6) .* B_cap(2:idx_c_6+1));
    atm_strike_6y = (1 - B_cap(idx_c_6+1)) / BPV_6y;

    BPV_10y = sum(delta_fwd(1:idx_c_10) .* B_cap(2:idx_c_10+1));
    atm_strike_10y = (1 - B_cap(idx_c_10+1)) / BPV_10y;

    % 3. Interpolate market flat_vols at the new ATM strikes
    % Obtain specific [n_maturities x 1] vectors for each Cap
    flat_atm_6y = zeros(length(maturities), 1);
    flat_atm_10y = zeros(length(maturities), 1);
    for m = 1:length(maturities)
        flat_atm_6y(m) = interp1(strikes, flat_vols(m, :), atm_strike_6y * 100, 'linear', 'extrap');
        flat_atm_10y(m) = interp1(strikes, flat_vols(m, :), atm_strike_10y * 100, 'linear', 'extrap');
    end

    % 4. Create the Vega Matrix (Central Difference)
    % Rows: Bucket 0-6y, Bucket 6-10y | Columns: Cap 6y, Cap 10y
    vega_caps = zeros(2, 2);
    bucket_rows = {1:6, 7:10};

    for b = 1:2
        rows = bucket_rows{b};

        % ==========================================
        % BUMP UP (+dVol)
        % ==========================================
        flat_up_6 = flat_atm_6y;   flat_up_6(rows)  = flat_up_6(rows) + dVol;
        flat_up_10 = flat_atm_10y; flat_up_10(rows) = flat_up_10(rows) + dVol;

        % Recalculate bumped spot vols using the LMM function (passing a single strike)
        [spot_up_6, ~, ~, ~, ~, ~, ~, ~] = lmm_spot_vols(flat_up_6, atm_strike_6y*100, maturities, dates, discounts, t0);
        [spot_up_10, ~, ~, ~, ~, ~, ~, ~] = lmm_spot_vols(flat_up_10, atm_strike_10y*100, maturities, dates, discounts, t0);

        % Robust substitution: loop through each individual caplet (UP)
        % UP price calculated by summing individual caplets (Black formula)
        price_up_6 = 0;
        for i = 1:idx_c_6
            price_up_6 = price_up_6 + caplet_black_LMM(fwd_rates(i), atm_strike_6y, delta_fwd(i), B_cap(i+1), tau_expiry(i), spot_up_6(i));
        end
        
        price_up_10 = 0;
        for i = 1:idx_c_10
            price_up_10 = price_up_10 + caplet_black_LMM(fwd_rates(i), atm_strike_10y, delta_fwd(i), B_cap(i+1), tau_expiry(i), spot_up_10(i));
        end

        % ==========================================
        % BUMP DOWN (-dVol)
        % ==========================================
        flat_dn_6 = flat_atm_6y;   flat_dn_6(rows)  = flat_dn_6(rows) - dVol;
        flat_dn_10 = flat_atm_10y; flat_dn_10(rows) = flat_dn_10(rows) - dVol;

        % Recalculate downward bumped spot vols
        [spot_dn_6, ~, ~, ~, ~, ~, ~, ~] = lmm_spot_vols(flat_dn_6, atm_strike_6y*100, maturities, dates, discounts, t0);
        [spot_dn_10, ~, ~, ~, ~, ~, ~, ~] = lmm_spot_vols(flat_dn_10, atm_strike_10y*100, maturities, dates, discounts, t0);

        % Robust substitution: loop through each individual caplet (DOWN)
        % DOWN price
        price_dn_6 = 0;
        for i = 1:idx_c_6
            price_dn_6 = price_dn_6 + caplet_black_LMM(fwd_rates(i), atm_strike_6y, delta_fwd(i), B_cap(i+1), tau_expiry(i), spot_dn_6(i));
        end
        
        price_dn_10 = 0;
        for i = 1:idx_c_10
            price_dn_10 = price_dn_10 + caplet_black_LMM(fwd_rates(i), atm_strike_10y, delta_fwd(i), B_cap(i+1), tau_expiry(i), spot_dn_10(i));
        end

        % ==========================================
        % VEGA CALCULATION (Matched to your formula)
        % ==========================================
        vega_caps(b, 1) = ((price_up_6 - price_dn_6) / 2) * N;
        vega_caps(b, 2) = ((price_up_10 - price_dn_10) / 2) * N;
    end
end