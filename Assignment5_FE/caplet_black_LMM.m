function price = caplet_black_LMM(L, K, delta, B_pay, tau, sigma)
% =========================================================
%  Black caplet formula under LMM (slide 20, Baviera)
%
%  caplet_i = B(0,T_{i+1}) * delta_i * [L_i*N(d1) - K*N(d2)]
%
%  with (slide 20):
%    d1 = (1/(sigma*sqrt(T_i))) * ln(L_i/K) + 0.5*sigma*sqrt(T_i)
%    d2 = d1 - sigma*sqrt(T_i)
%
%  INPUTS:
%    L      : forward rate L_i(t0)
%    K      : strike
%    delta  : Act/360 year fraction of the period
%    B_pay  : discount factor to payment date B(0, T_{i+1})
%    tau    : Act/365 time to expiry (= T_i - t0)/365
%    sigma  : spot vol of this caplet
%
%  OUTPUT:
%    price  : caplet price (per unit notional)
% =========================================================

    % Handle degenerate cases
    if sigma <= 0 || tau <= 0 || L <= 0 || K <= 0
        price = max(delta * B_pay * (L - K), 0);
        return;
    end
    
    vol_sqrt_T = sigma * sqrt(tau);
    
    % d1 and d2 from Black formula (slide 20, LMM)
    d1 = (log(L/K) + 0.5 * sigma^2 * tau) / vol_sqrt_T;
    d2 = d1 - vol_sqrt_T;
    
    % Black caplet price (slide 6 and 20)
    price = B_pay * delta * (L * normcdf(d1) - K * normcdf(d2));
end