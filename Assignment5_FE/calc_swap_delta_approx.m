function delta_swap = calc_swap_delta_approx(t0, T_expiry, dates_df, df_base, df_bump)
% CALC_SWAP_DELTA_APPROX Computes Swap PV01 using the Trader's Approximation.
% Assumes the floating leg is perfectly at par (Delta = 0) and all interest 
% rate risk comes solely from the fixed leg (Annuity/BPV).

    % 1. Auto-generate payment dates (Modified Following)
    maturity_years = round(yearfrac(t0, T_expiry, 3));
    dates_swap_raw = datemnth(t0, 12 * (1:maturity_years)');
    dates_swap = busdate(dates_swap_raw, 'modifiedfollow');

    % 2. Calculate year fractions for Fixed Leg (30E/360 basis = 6)
    dates_all = [t0; dates_swap(1:end-1)];
    delta_t = yearfrac(dates_all, dates_swap, 6);

    % 3. Interpolate DFs at payment dates
    df_swap_base = linearRateInterp(dates_df, df_base, t0, dates_swap);
    df_swap_bump = linearRateInterp(dates_df, df_bump, t0, dates_swap);

    % ==========================================
    % TRADER'S APPROXIMATION LOGIC
    % ==========================================
    
    % A. Calcolo il BPV e il Tasso Par (S) sulla curva base
    BPV_base = sum(delta_t .* df_swap_base);
    float_leg_base = 1 - df_swap_base(end);
    S = float_leg_base / BPV_base;
    
    % B. Calcolo solo il BPV sulla curva bumpata
    BPV_bump = sum(delta_t .* df_swap_bump);
    
    % C. Calcolo del Delta (Posizione Payer)
    % Ignoriamo la gamba variabile. 
    % Se i tassi salgono, il BPV scende. Il valore attuale di ciò che DEVO pagare 
    % (S * BPV) diminuisce, generando un profitto (Delta positivo).
    delta_swap = (S * BPV_base) - (S * BPV_bump);

end