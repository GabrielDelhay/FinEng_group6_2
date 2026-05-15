function X = price_certificate3y(p1, p12, p3, N, c1, c2, B1, B2, B3, BPV_y1, BPV_y2, BPV_y3, spread)
% Prices a 3y autocallable certificate, extending price_certificate to a third
% possible redemption date T3. Early redemption triggers at T1 or T2 if S < K.
%
% INPUTS:
    % p1     : P(ER at T1) = P(S(T1) < K)
    % p12    : P(survive T1, ER at T2) = P(S(T1) >= K, S(T2) < K)
    % p3     : P(reach T3) = P(S(T1) >= K, S(T2) >= K)
    % N      : notional
    % c1     : coupon paid at T1 or T2 if early redemption triggers
    % c2     : coupon paid at T3 if no early redemption
    % B1, B2, B3   : discount factors B(t0, T1), B(t0, T2), B(t0, T3)
    % BPV_y1, BPV_y2, BPV_y3 : BPV of floating leg in each year (Act/360, quarterly)
    % spread : Euribor3M spread on the floating leg
%
% OUTPUT:
    % X      : upfront payment as a fraction of notional
    NPV_float   = N*(1-B1) + spread*N*BPV_y1 ...
                + (1-p1) * (N*(B1-B2) + spread*N*BPV_y2) ...
                + p3     * (N*(B2-B3) + spread*N*BPV_y3);
    NPV_coupons = N * (c1*B1*p1 + c1*B2*p12 + c2*B3*p3);
    X           = (NPV_float - NPV_coupons) / N;
end