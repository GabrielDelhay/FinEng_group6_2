% Price caplets in segment [i_alpha, i_beta] with linearly interpolated BMM spot vol.
% sa, sb : left and right vol boundaries for the segment
% T_alpha, T_beta : reset dates at segment boundaries (used for linear interp weight)
function p = priceBMM_segment(sb, sa, i_alpha, i_beta, B0_T, delta, L0, ...
                               resetDates, t0, T_alpha, T_beta)
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