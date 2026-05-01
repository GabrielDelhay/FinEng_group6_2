function price = capletBMM_strike(nu, Bi1, di, Li, K, Ti_t0)
% BMM caplet price at arbitrary strike K
    if nu * sqrt(Ti_t0) < 1e-10
        price = 0; return
    end
    d1    = log((1 + di*Li)/(1 + di*K)) / (nu*sqrt(Ti_t0)) + 0.5*nu*sqrt(Ti_t0);
    d2    = d1 - nu * sqrt(Ti_t0);
    price = Bi1 * ((1 + di*Li)*normcdf(d1) - (1 + di*K)*normcdf(d2));
end