function V = ZBPutHW(dates, discounts, t0, T_exp, T_mat, Kstrike, a, sigma)
% ZBPutHW  Closed-form put on a zero-coupon bond in 1-factor Hull-White.
%
% Payoff at T_exp:  max( Kstrike - P(T_exp, T_mat),  0 ).
%
% Pricing formula (Black-style under T_exp-forward measure):
%
%     V = Kstrike * B0(T_exp) * Phi(-h + sigma_p)  -  B0(T_mat) * Phi(-h)
%
%     sigma_p = sigma * sqrt((1 - exp(-2*a*T_exp))/(2*a)) * Bhw(T_exp, T_mat)
%     h       = (1/sigma_p) * log( B0(T_mat) / (B0(T_exp) * Kstrike) ) + sigma_p/2
%
% NOTATION
%   B0(.) : initial discount factor (ZCB at t=0) from the market curve.
%
% VECTORIZATION
%   T_exp   : scalar expiry date.
%   T_mat   : scalar or column vector of bond maturity dates.
%   Kstrike : scalar or column vector, same size as T_mat.
%   The output V has the same size as T_mat.
%
% This helper is the shared building block of both capletHW and
% europeanPayerSwaptionHW (Jamshidian decomposition).

    %% Year fractions
    tau_exp = yearfrac(t0, T_exp, 3);
    tau_mat = yearfrac(t0, T_mat, 3);
    tau_mat = tau_mat(:);
    Kstrike = Kstrike(:);

    %% Initial discount factors
    B0_exp = linearRateInterp(dates, discounts, t0, T_exp);
    B0_mat = arrayfun(@(d) linearRateInterp(dates, discounts, t0, d), T_mat);
    B0_mat = B0_mat(:);

    %% HW bond volatility seen from t=0
    Bhw     = (1 - exp(-a*(tau_mat - tau_exp))) / a;
    sigma_p = sigma * sqrt((1 - exp(-2*a*tau_exp))/(2*a)) .* Bhw;

    %% Black-style put on the ZCB
    h = (1 ./ sigma_p) .* log(B0_mat ./ (B0_exp .* Kstrike)) + sigma_p/2;
    V = Kstrike .* B0_exp .* normcdf(-h + sigma_p) - B0_mat .* normcdf(-h);
end
