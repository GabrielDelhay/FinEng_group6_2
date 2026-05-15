function p = digital_put_NIG(p_opt, alpha, x_K, tau)
% Gil-Pelaez inversion
% INPUTS:
    % p_opt  : calibrated NIG parameters [theta, kappa, sigma]
    % alpha  : mixing exponent = 0.5 (NIG)
    % x_K    : log-moneyness log(K/F0) at reset date
    % tau    : time to reset date in years
%
% OUTPUT:
    % p      : risk-neutral probability P(S(tau) < K) via Gil-Pelaez inversion:
    %          p = 1/2 - (1/pi) * int_0^inf Im[e^{-i*u*x} * phi(u)] / u du
    phi = @(u) char_fun(u, p_opt, alpha, tau);
    p   = 0.5 - (1/pi) * integral(@(u) imag(exp(-1i*u*x_K).*phi(u))./u, ...
                                   1e-8, 500, 'AbsTol',1e-10,'RelTol',1e-8);
end
