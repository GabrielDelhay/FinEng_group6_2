function price = capletBMM(nu, Bi1, di, Li, Ti_t0)
% BMM caplet price at ATM strike K = Li (slide 20, Baviera 2006).
    if nu * sqrt(Ti_t0) < 1e-10
        price = 0; return
    end
    d1    =  0.5 * nu * sqrt(Ti_t0);
    d2    = -0.5 * nu * sqrt(Ti_t0);
    price = Bi1 * (1 + di*Li) * (normcdf(d1) - normcdf(d2));
end