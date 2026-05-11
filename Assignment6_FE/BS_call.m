function C = BS_call(F0, K, B, sigma, T)
    % Black (1976) call price formula 
    % c(K) = B(t0,t) * { F0*N[d1] - K*N[d2] }
    % with d_{1,2} = ln(F0/K)/sqrt(T*sigma^2) +/- (1/2)*sqrt(T*sigma^2)
    d1 = (log(F0./K) + 0.5*sigma.^2.*T) ./ (sigma.*sqrt(T));
    d2 = d1 - sigma.*sqrt(T);
    C  = B .* (F0 .* normcdf(d1) - K .* normcdf(d2));
end