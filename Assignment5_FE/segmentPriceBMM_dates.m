function price = segmentPriceBMM_dates(sigma_beta, sigma_alpha, ...
                  i_alpha, i_beta, discGrid, delta, L0, resetDates, t0, ...
                  T_alpha, T_beta)
% Price caplets in [i_alpha, i_beta] with linearly interpolated vols
% using reset dates as interpolation variable (consistent with market practice)
    price = 0;
    for i = i_alpha:i_beta
        w     = (resetDates(i+1) - T_alpha) / (T_beta - T_alpha);
        nu_i  = sigma_alpha + w * (sigma_beta - sigma_alpha);
        Ti_t0 = yearfrac(t0, resetDates(i+1), 3);
        price = price + capletBMM_strike(nu_i, discGrid(i+2), delta(i+1), ...
                                          L0(i+1), L0(i+1), Ti_t0);
    end
end