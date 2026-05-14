function V = capletHW(dates, discounts, t0, T_reset, T_pay, K, a, sigma)
% Closed-form caplet in 1-factor Hull-White.
%
% Payoff at T_pay:  delta * max( L(T_reset, T_pay) - K, 0 ).
%
% Identity used:
%     Caplet  =  (1 + K*delta) * Put on ZCB with strike 1/(1 + K*delta)
% The ZCB put has a Black-style closed form, implemented by ZBPutHW.
%
% NOTATION
%   B0(.) : initial discount factor (ZCB at t=0), from market curve.

    %% Accrual fraction
    tau_reset = yearfrac(t0, T_reset, 3);
    tau_pay   = yearfrac(t0, T_pay,   3);
    delta     = tau_pay - tau_reset;

    %% Strike of the underlying ZB put
    K_zbp = 1 / (1 + K*delta);

    %% Caplet = (1 + K*delta) * ZBPut at strike K_zbp
    V = (1 + K*delta) * ZBPutHW(dates, discounts, t0, T_reset, T_pay, K_zbp, a, sigma);
end
