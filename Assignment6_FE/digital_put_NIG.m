%% digital_put_NIG.m  — Gil-Pelaez inversion
function p = digital_put_NIG(p_opt, alpha, x_K, tau)
    phi = @(u) char_fun(u, p_opt, alpha, tau);
    p   = 0.5 - (1/pi) * integral(@(u) imag(exp(-1i*u*x_K).*phi(u))./u, ...
                                   1e-8, 500, 'AbsTol',1e-10,'RelTol',1e-8);
end