%% digital_put_Black.m
function p = digital_put_Black(F, K, sigma, tau)
    d2 = (log(F/K) - 0.5*sigma^2*tau) / (sigma*sqrt(tau));
    p  = normcdf(-d2);
end