function bucket_delta = coarse_bucket_delta(ratesSet, datesSet, dates, flat_vols, strikes, maturities, t0, N, spread, bond, true_price, BPV, n_depos, n_futures)
bucket_delta = zeros(3, 1);

% --- Bucket 1: 0-2y ---
% Bump all depos and futures
ratesSet.depos(1:n_depos, :)     = ratesSet.depos(1:n_depos, :)     + BPV;
ratesSet.futures(1:n_futures, :) = ratesSet.futures(1:n_futures, :) + BPV;
% Bump swap_2 (index 2)
ratesSet.swaps(2, :) = ratesSet.swaps(2, :) + BPV;
[~, disc_bump, ~] = bootstrap(datesSet, ratesSet);
[spot_bump, ~, B_bump, fwd_bump, delta_bump, tau_bump, ~, idx_bump] = ...
    lmm_spot_vols(flat_vols, strikes, maturities, dates, disc_bump, t0);
[~, ~, X_flat_bump, ~, ~] = price_structured_bond(N, spread, bond, B_bump, delta_bump, tau_bump, fwd_bump, idx_bump, spot_bump, strikes);
bucket_delta(1) = X_flat_bump * N - true_price;

% Restore
ratesSet.depos(1:n_depos, :)     = ratesSet.depos(1:n_depos, :)     - BPV;
ratesSet.futures(1:n_futures, :) = ratesSet.futures(1:n_futures, :) - BPV;
ratesSet.swaps(2, :)             = ratesSet.swaps(2, :)             - BPV;

% --- Bucket 2: 2-6y ---
% Bump swaps idx 3 to 6
ratesSet.swaps(3:6, :) = ratesSet.swaps(3:6, :) + BPV;
[~, disc_bump, ~] = bootstrap(datesSet, ratesSet);
[spot_bump, ~, B_bump, fwd_bump, delta_bump, tau_bump, ~, idx_bump] = ...
    lmm_spot_vols(flat_vols, strikes, maturities, dates, disc_bump, t0);
[~, ~, X_flat_bump, ~, ~] = price_structured_bond(N, spread, bond, B_bump, delta_bump, tau_bump, fwd_bump, idx_bump, spot_bump, strikes);
bucket_delta(2) = X_flat_bump * N - true_price;

% Restore
ratesSet.swaps(3:6, :) = ratesSet.swaps(3:6, :) - BPV;

% --- Bucket 3: 6-10y ---
% Bump swaps idx 7 to 10
ratesSet.swaps(7:10, :) = ratesSet.swaps(7:10, :) + BPV;
[~, disc_bump, ~] = bootstrap(datesSet, ratesSet);
[spot_bump, ~, B_bump, fwd_bump, delta_bump, tau_bump, ~, idx_bump] = ...
    lmm_spot_vols(flat_vols, strikes, maturities, dates, disc_bump, t0);
[~, ~, X_flat_bump, ~, ~] = price_structured_bond(N, spread, bond, B_bump, delta_bump, tau_bump, fwd_bump, idx_bump, spot_bump, strikes);
bucket_delta(3) = X_flat_bump * N - true_price;

% Restore
ratesSet.swaps(7:10, :) = ratesSet.swaps(7:10, :) - BPV;

end