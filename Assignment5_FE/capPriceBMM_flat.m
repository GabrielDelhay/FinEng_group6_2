function price = capPriceBMM_flat(flatVol, discGrid, delta, L0, resetDates, t0, nCaplets)
% Price of a cap with nCaplets caplets, all priced with the same flat vol.
    price = 0;
    for i = 1:nCaplets
        Ti_t0 = yearfrac(t0, resetDates(i+1), 3);
        Bi1   = discGrid(i+2);
        di    = delta(i+1);
        Li    = L0(i+1);
        price = price + capletBMM(flatVol, Bi1, di, Li, Ti_t0);
    end
end