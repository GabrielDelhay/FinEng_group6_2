function bucket_vega = coarse_bucket_vega(flat_vols, strikes, maturities, dates, discounts, ...
                                           t0, N, spread, bond, B_cap, delta_fwd, tau_expiry, ...
                                           fwd_rates, cap_maturity_idx, dVol, cap_dates)
% coarse_bucket_vega  Bucketed vega (0-6y and 6-10y) by bumping flat_vols columns.
%
%   Bucket 1: 0-6y  -> maturities 1,2,3,4,5,6  -> columns 1:6
%   Bucket 2: 6-10y -> maturities 7,8,9,10      -> columns 7:10
%
%   Uses central finite difference: (price_up - price_down) / (2 * dVol)
%   scaled to a 1% vol move.

bucket_vega = zeros(2, 1);

% Index ranges in the maturities vector
% maturities = [1,2,3,4,5,6,7,8,9,10,12,15,20]
% 0-6y  -> indices 1:6  (maturities 1y to 6y)
% 6-10y -> indices 7:10 (maturities 7y to 10y)
bucket_rows = {1:6, 7:10};

for b = 1:2
    rows = bucket_rows{b};

    % --- Up bump ---
    flat_vols_up = flat_vols;
    flat_vols_up(rows,:) = flat_vols_up(rows,:) + dVol;
    [spot_vols_up, ~, ~, ~, ~, ~, ~, ~] = ...
        lmm_spot_vols(flat_vols_up, strikes, maturities, dates, discounts, t0);
    [~, ~, X_up, ~, ~] = price_structured_bond(N, spread, bond, B_cap, delta_fwd, tau_expiry, fwd_rates, cap_maturity_idx, spot_vols_up, strikes, cap_dates, t0);

    % --- Down bump ---
    flat_vols_down = flat_vols;
    flat_vols_down(rows,:) = flat_vols_down(rows,:) - dVol;
    [spot_vols_down, ~, ~, ~, ~, ~, ~, ~] = ...
        lmm_spot_vols(flat_vols_down, strikes, maturities, dates, discounts, t0);
    [~, ~, X_down, ~, ~] = price_structured_bond(N, spread, bond, B_cap, delta_fwd, tau_expiry, ...
                                                  fwd_rates, cap_maturity_idx, spot_vols_down, strikes, cap_dates, t0);

    % Central difference, scaled to +1% vol move (same convention as total vega)
    bucket_vega(b) = (X_up - X_down) * N / 2;
end
end