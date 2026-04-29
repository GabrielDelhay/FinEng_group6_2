function price = segmentPrice(nu_right, nu_left, i_left, i_right, ...
                               discGrid, delta, L0, resetDates, t0)
% Price of caplets in [i_left, i_right] with linearly interpolated spot vols.
    price = 0;
    for i = i_left:i_right
        nu_i  = nu_left + (i - i_left) / (i_right - i_left + 1) * (nu_right - nu_left);
        Ti_t0 = yearfrac(t0, resetDates(i+1), 3);
        Bi1   = discGrid(i+2);
        di    = delta(i+1);
        Li    = L0(i+1);
        price = price + capletBMM(nu_i, Bi1, di, Li, Ti_t0);
    end
end