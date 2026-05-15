function X = price_certificate(p1, N, c1, c2, B1, B2, BPV_y1, BPV_y2, spread)
% Prices a 2y autocallable certificate by computing the upfront X as a % of notional.
% The structure redeems early at T1 if S(T1_reset) < K (probability p1),
% otherwise runs to T2. Funding leg is Euribor3M + spread, coupons are c1 (if ER) or c2.
%
% INPUTS:
    % p1      : risk-neutral probability of early redemption at T1, P(S(T1) < K)
    % N       : notional
    % c1      : coupon paid at T1 if early redemption triggers
    % c2      : coupon paid at T2 if no early redemption
    % B1, B2  : discount factors B(t0, T1) and B(t0, T2)
    % BPV_y1  : BPV of floating leg in year 1 (Act/360, quarterly)
    % BPV_y2  : BPV of floating leg in year 2 (Act/360, quarterly)
    % spread  : Euribor3M spread on the floating leg
% OUTPUT:
    % X       : upfront payment as a fraction of notional (NPV_float - NPV_coupons) / N
    NPV_float   = N*(1-B1) + spread*N*BPV_y1 ...
                + (1-p1)*(N*(B1-B2) + spread*N*BPV_y2);
    NPV_coupons = N*(c1*B1*p1 + c2*B2*(1-p1));
    X           = (NPV_float - NPV_coupons) / N;
end