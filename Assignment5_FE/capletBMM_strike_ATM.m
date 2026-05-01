% BMM ATM caplet price. At-the-money (K = Li1) exact simplification.
function price = capletBMM_strike_ATM(nu, Bi2, di1, Li1, Ti_t0)
    % nu    : BMM spot vol for this caplet (vol of log forward bond)
    % Bi2   : B0_T(i+2) = P(t0, T_{i+1}), discount to payment date T_{i+1}
    % di1   : delta(i+1), year fraction of the accrual period [T_i, T_{i+1}]
    % Li1   : L0(i+1), ATM forward Libor = strike
    % Ti_t0 : yearfrac(t0, T_{i+1}), time used for vol scaling (consistent with bootstrap)
    if nu * sqrt(Ti_t0) < 1e-10
        price = 0; return
    end
    % ATM exact: d1 = +0.5*nu*sqrt(T), d2 = -d1
    % Price = Bi2*(1+di1*Li1)*[N(d1)-N(d2)] = Bi2*(1+di1*Li1)*[2*N(d1)-1]
    d1    = 0.5 * nu * sqrt(Ti_t0);
    price = Bi2 * (1 + di1*Li1) * (2*normcdf(d1) - 1);
end