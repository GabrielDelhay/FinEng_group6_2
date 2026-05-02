function flatVolsATM = getATMFlatVols(B0_T, delta, maturitiesYears, strikes_table, vols_table)
    K_ATM = zeros(1, numel(maturitiesYears));
    for m = 1:numel(maturitiesYears)
        nQ       = 4 * m;
        K_ATM(m) = (1 - B0_T(nQ+1)) / sum(delta(1:nQ) .* B0_T(2:nQ+1));
    end
    flatVolsATM = zeros(1, numel(maturitiesYears));
    for m = 1:numel(maturitiesYears)
        flatVolsATM(m) = interp1(strikes_table, vols_table(m,:), K_ATM(m), 'linear', 'extrap');
    end
end
