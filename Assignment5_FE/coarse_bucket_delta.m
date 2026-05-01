function bucket_delta = calc_macro_buckets(delta, n_depos, n_futures)
% CALC_MACRO_BUCKETS Aggregates pointwise Deltas into 3 macro-buckets (0-2y, 2-6y, 6-10y)
% using a simple linear sum.
%
% Inputs:
%   delta     : Vector of pointwise PV01s (e.g., length 24) from Ex 1.c
%   n_depos   : Number of bumped deposits (e.g., 3)
%   n_futures : Number of bumped futures (e.g., 7)

    bucket_delta = zeros(3, 1);
    
    % Index where depos and futures end, and swaps begin
    idx_swaps_start = n_depos + n_futures; 
    
    % --- Bucket 1: 0-2y ---
    % Includes: all depos (1-3), all futures (1-7), and 2y Swap (1st saved swap)
    idx_b1 = 1 : (idx_swaps_start + 1);
    bucket_delta(1) = sum(delta(idx_b1));
    
    % --- Bucket 2: 2-6y ---
    % Includes: Swaps 3, 4, 5, 6 (2nd to 5th swap in the vector)
    idx_b2 = (idx_swaps_start + 2) : (idx_swaps_start + 5);
    bucket_delta(2) = sum(delta(idx_b2));
    
    % --- Bucket 3: 6-10y ---
    % Includes: Swaps 7, 8, 9, 10 (6th to 9th swap in the vector)
    % Note: Stopping at 10y aligns with the bond maturity.
    idx_b3 = (idx_swaps_start + 6) : (idx_swaps_start + 9);
    bucket_delta(3) = sum(delta(idx_b3));

end