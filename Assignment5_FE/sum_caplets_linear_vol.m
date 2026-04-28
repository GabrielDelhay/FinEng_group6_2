function total_price = sum_caplets_linear_vol(i_a, i_b, ...
    fwd_rates, K, delta_fwd, all_B, tau_expiry, ...
    sigma_alpha, sigma_beta, T_alpha, T_beta)

% =========================================================
%  Sum of caplet prices within a bucket, where spot vols
%  are LINEARLY interpolated between sigma_alpha and sigma_beta
%  (slide 23, linear constraint)
%
%  INPUTS:
%    i_a, i_b       : first and last caplet index in bucket
%    fwd_rates      : vector of forward rates
%    K              : strike
%    delta_fwd      : vector of Act/360 year fractions
%    all_B          : vector of discount factors (indexed i+1 for payment)
%    tau_expiry     : vector of Act/365 times to expiry
%    sigma_alpha    : spot vol at left edge of bucket (known)
%    sigma_beta     : spot vol at right edge of bucket (unknown)
%    T_alpha, T_beta: tau_expiry values at bucket edges
%
%  OUTPUT:
%    total_price : sum of caplet prices in this bucket
% =========================================================

    total_price = 0;
    
    for i = i_a:i_b
        % Linear interpolation of spot vol (slide 23)
        if T_beta > T_alpha
            sigma_i = sigma_alpha + ...
                (tau_expiry(i) - T_alpha) / (T_beta - T_alpha) ...
                * (sigma_beta - sigma_alpha);
        else
            sigma_i = sigma_beta; % single caplet, no interpolation needed
        end
        
        % Ensure vol is positive (fzero may try negative values)
        sigma_i = max(sigma_i, 1e-8);
        
        total_price = total_price + ...
            caplet_black_LMM(fwd_rates(i), K, delta_fwd(i), ...
                              all_B(i+1), tau_expiry(i), sigma_i);
    end
end