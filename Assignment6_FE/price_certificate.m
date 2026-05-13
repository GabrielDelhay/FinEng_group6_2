%% price_certificate.m
function X = price_certificate(p1, N, c1, c2, B1, B2, BPV_y1, BPV_y2, spread)
    NPV_float   = N*(1-B1) + spread*N*BPV_y1 ...
                + (1-p1)*(N*(B1-B2) + spread*N*BPV_y2);
    NPV_coupons = N*(c1*B1*p1 + c2*B2*(1-p1));
    X           = (NPV_float - NPV_coupons) / N;
end