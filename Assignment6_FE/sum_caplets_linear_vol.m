function total_price = sum_caplets_linear_vol(i_a, i_b, ...
    fwd_rates, K, delta_fwd, all_B, T_reset, ...
    sigma_alpha, sigma_beta, T_alpha, T_beta, tau)

% SUM_CAPLETS_LINEAR_VOL  Total price of the caplets in a bucket, with
% spot vols linearly interpolated between the left edge sigma_alpha
% (known) and the right edge sigma_beta (unknown).
%
%  INPUTS:
%    i_a, i_b      : first and last caplet index in the bucket
%    fwd_rates     : vector of forward rates per caplet
%    K             : strike
%    delta_fwd     : vector of Act/360 year fractions per caplet
%    all_B         : discount factors at caplet payment dates
%    T_reset       : reset times of every caplet (used for interpolation)
%    sigma_alpha   : spot vol at the left edge of the bucket (known)
%    sigma_beta    : spot vol at the right edge of the bucket (unknown)
%    T_alpha       : T_reset at the left edge of the bucket
%    T_beta        : T_reset at the right edge of the bucket
%    tau           : Act/365 times to caplet reset (Black scaling)
%
%  OUTPUT:
%    total_price : sum of caplet prices in this bucket

    total_price = 0;
    
    for i = i_a:i_b
        % Linear interpolation of the spot vol at this caplet.

            sigma_i = sigma_alpha + ...
                (T_reset(i) - T_alpha) / (T_beta - T_alpha) ...
                * (sigma_beta - sigma_alpha);
        
            total_price = total_price + ...
                caplet_black_LMM(fwd_rates(i), K, delta_fwd(i), ...
                                  all_B(i+1), tau(i), sigma_i);
    end

end