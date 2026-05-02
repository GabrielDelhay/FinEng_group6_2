function [delta, T_delta] = calc_delta_buckets(N, spread, bond, datesSet, ratesSet, flat_vols, strikes, maturities, t0, true_price,n_depos, n_futures, n_swaps, n_total, BPV, cap_dates)

    price = zeros(n_total, 1);
    idx = 1; 

    % 1. Bump Depos
    for i = 1:n_depos
        ratesSet.depos(i, :) = ratesSet.depos(i, :) + BPV;  % Bump
        [d_bump, disc_bump, ~] = bootstrap(datesSet, ratesSet);
        
        [spot_bump, ~, B_bump, fwd_bump, delta_bump, tau_bump, ~, idx_bump] = ...
            lmm_spot_vols(flat_vols, strikes, maturities, d_bump, disc_bump, t0);
        
        % Extract X_flat (3rd output)
        [~, ~, X_flat_bump, ~, ~] = price_structured_bond(N, spread, bond, B_bump, delta_bump, tau_bump, fwd_bump, idx_bump, spot_bump, strikes, cap_dates, t0);
        price(idx) = X_flat_bump * N;
        
        ratesSet.depos(i, :) = ratesSet.depos(i, :) - BPV;  % Restore
        idx = idx + 1;
    end

    % 2. Bump Futures
    for i = 1:n_futures
        ratesSet.futures(i, :) = ratesSet.futures(i, :) + BPV; 
        [d_bump, disc_bump, ~] = bootstrap(datesSet, ratesSet);
        
        [spot_bump, ~, B_bump, fwd_bump, delta_bump, tau_bump, ~, idx_bump] = ...
            lmm_spot_vols(flat_vols, strikes, maturities, d_bump, disc_bump, t0);
            
        [~, ~, X_flat_bump, ~, ~] = price_structured_bond(N, spread, bond, B_bump, delta_bump, tau_bump, fwd_bump, idx_bump, spot_bump, strikes, cap_dates, t0);
        price(idx) = X_flat_bump * N;
        
        ratesSet.futures(i, :) = ratesSet.futures(i, :) - BPV;
        idx = idx + 1;
    end

    % 3. Bump Swaps (skip the first)
    for i = 2:n_swaps
        ratesSet.swaps(i, :) = ratesSet.swaps(i, :) + BPV; 
        [d_bump, disc_bump, ~] = bootstrap(datesSet, ratesSet);
        
        [spot_bump, ~, B_bump, fwd_bump, delta_bump, tau_bump, ~, idx_bump] = ...
            lmm_spot_vols(flat_vols, strikes, maturities, d_bump, disc_bump, t0);
            
        [~, ~, X_flat_bump, ~, ~] = price_structured_bond(N, spread, bond, B_bump, delta_bump, tau_bump, fwd_bump, idx_bump, spot_bump, strikes, cap_dates, t0);
        price(idx) = X_flat_bump * N;
        
        ratesSet.swaps(i, :) = ratesSet.swaps(i, :) - BPV; 
        idx = idx + 1;
    end

    % Vectorized delta calculation
    delta = price - true_price;

    % Create labels and table
    labels = [strcat("Depo_", string(1:n_depos)), ...
              strcat("Future_", string(1:n_futures)), ...
              strcat("Swap_", string(2:n_swaps))].';
              
    T_delta = table(labels, delta, 'VariableNames', {'Pillar', 'PV01_EUR'});
end