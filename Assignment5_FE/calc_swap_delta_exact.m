function delta_swap = calc_swap_delta_exact(t0, T_expiry, dates_df, df_base, df_bump)
% CALC_SWAP_DELTA_EXACT Computes Swap PV01 generating payment dates internally.
%
% Inputs:
%   t0         : Settlement date (swap start date)
%   T_expiry   : Expiry date of the swap
%   dates_df   : Discount factor dates from the bootstrap curve
%   df_base    : Baseline discount factors
%   df_bump    : Bumped discount factors

    % 1. Auto-generate payment dates
    % Calculate maturity in years (rounded to nearest integer), convention
    % doesn't come into play
    maturity_years = round(yearfrac(t0, T_expiry, 3));
    
    % Generate annual payment dates with Modified Following convention
    dates_swap_raw = datemnth(t0, 12 * (1:maturity_years)');
    dates_swap = busdate(dates_swap_raw, 'modifiedfollow');

    % 2. Calculate year fractions (30E/360 basis = 6)
    dates_all = [t0; dates_swap(1:end-1)];
    delta_t = yearfrac(dates_all, dates_swap, 6);

    % 3. Interpolate DFs at payment dates
    df_swap_base = linearRateInterp(dates_df, df_base, t0, dates_swap);
    df_swap_bump = linearRateInterp(dates_df, df_bump, t0, dates_swap);

    % ==========================================
    % BASELINE SCENARIO (Par Rate)
    % ==========================================
    BPV_base = sum(delta_t .* df_swap_base);
    float_leg_base = 1 - df_swap_base(end);
    S = float_leg_base / BPV_base;

    % ==========================================
    % BUMPED SCENARIO (Revaluation)
    % ==========================================
    BPV_bump = sum(delta_t .* df_swap_bump);
    float_leg_bump = 1 - df_swap_bump(end);
    
    NPV_fixed_bump = S * BPV_bump;
    NPV_float_bump = float_leg_bump;
    
    % ==========================================
    % DELTA CALCULATION (Payer Position)
    % ==========================================
    delta_swap = NPV_float_bump - NPV_fixed_bump;

end