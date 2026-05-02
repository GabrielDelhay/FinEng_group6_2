function p = priceBMM_segment(sb, sa, i_alpha, i_beta, B0_T, delta, L0, ...
                               resetDates, t0, T_alpha, T_beta)
% Prices the caplets in index range [i_alpha, i_beta] using a BMM spot vol
% linearly interpolated between sa (left boundary) and sb (right boundary).
% Inputs: sb, sa (right/left vol boundaries), 
%         i_alpha, i_beta (caplet indices),
%         B0_T, 
%         delta, 
%         L0, 
%         resetDates, 
%         t0, 
%         T_alpha, T_beta (boundary dates)

    p = 0;
    for i = i_alpha:i_beta
        % Linear interpolation weight in [0,1] across the segment
        w     = (resetDates(i+1) - T_alpha) / (T_beta - T_alpha);
        nu_i  = sa + w * (sb - sa);
        Ti_t0 = yearfrac(t0, resetDates(i+2), 3);
        p     = p + capletBMM_strike_ATM(nu_i, B0_T(i+2), delta(i+1), ...
                                          L0(i+1), Ti_t0);
    end
end