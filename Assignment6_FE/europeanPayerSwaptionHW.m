function V = europeanPayerSwaptionHW(dates, discounts, t0, T_alpha, T_pay, K, a, sigma)
% European payer swaption - Jamshidian closed form, 1-factor Hull-White.
%
% Jamshidian's trick decomposes the swaption into a weighted sum of puts
% on individual zero-coupon bonds, exploiting the comonotonicity of all
% P(T_alpha, T_k; x) in the single state variable x.
%
% NOTATION
%   B0(T)        : initial ZCB price (discount factor) observed at t=0
%   P0(T_1..T_n) : initial coupon bond price = sum_k c_k * B0(T_k)
%   x*           : critical state s.t. the forward coupon bond at T_alpha = 1
%   K_k          : ZCB strike of the k-th put = P(T_alpha, T_k; x*)
%
% The actual put-on-ZCB pricing is delegated to ZBPutHW.

    %% Year fractions
    tau_alpha = yearfrac(t0, T_alpha, 3);
    tau_pay   = yearfrac(t0, T_pay,   3);
    tau_pay   = tau_pay(:);                                  % column

    %% Coupon coefficients c_k (notional = 1)
    deltas = [tau_pay(1) - tau_alpha; diff(tau_pay)];
    c      = K * deltas;
    c(end) = c(end) + 1;                                     % last cf includes notional

    %% Initial discount factors (ZCB at t=0)
    B0_alpha = linearRateInterp(dates, discounts, t0, T_alpha);
    B0_pay   = arrayfun(@(d) linearRateInterp(dates, discounts, t0, d), T_pay);
    B0_pay   = B0_pay(:);

    %% Initial coupon bond P0 (sanity-check quantity, equals 1 at x*)
    P0_fixed = sum(c .* B0_pay);                             %#ok<NASGU>

    %% A(T_alpha, T_k) and Bhw(T_alpha, T_k) - coefficients of the HW bond formula
    Bhw_k = (1 - exp(-a*(tau_pay - tau_alpha))) / a;
    term1 = (2*(1 - exp(-a*tau_alpha))/a)    .* (1 - exp(-a*(tau_pay - tau_alpha)));
    term2 = ((1 - exp(-2*a*tau_alpha))/(2*a)).* (1 - exp(-2*a*(tau_pay - tau_alpha)));
    V_k   = (sigma^2/a^2) * (term1 - term2);
    A_k   = (B0_pay / B0_alpha) .* exp(-0.5 * V_k);

    %% Find x* such that the forward coupon bond at T_alpha equals 1
    objective = @(x) sum(c .* A_k .* exp(-Bhw_k * x)) - 1;
    x_star    = fzero(objective, [-1, 1]);

    %% Strikes of the individual ZB puts:  K_k = P(T_alpha, T_k; x*)
    K_k = A_k .* exp(-Bhw_k * x_star);

    %% Sum of weighted ZB puts (closed form delegated to ZBPutHW)
    ZBPut = ZBPutHW(dates, discounts, t0, T_alpha, T_pay, K_k, a, sigma);
    V     = sum(c .* ZBPut);
end
