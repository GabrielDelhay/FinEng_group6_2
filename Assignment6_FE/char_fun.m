function phi = char_fun(xi, p, alpha, T)
    % Characteristic function of the log-return ft = ln(Ft/F0) 
    %
    % phi(xi) = exp{ -i*xi * ln L[eta] } * L[ (xi^2 + i*(1+2*eta)*xi) / 2 ]
    %
    % where:
    %   - exp{-i*xi * ln L[eta]} is the NORMALIZATION term 
    %     it ensures E[e^{ft}] = 1 (martingale property)
    %   - L[(xi^2 + i*(1+2*eta)*xi)/2] is the main Laplace transform
    %     evaluated at a COMPLEX argument (this is why we need analyticity)
    %
    % Note: ln L[eta] is real since eta is real
    
    eta = p(3);
    
    % Normalization: ln L evaluated at real point eta
    lnL_eta = laplace_exponent(eta, p, alpha, T);
    
    % Main argument: complex number (xi^2 + i*(1+2*eta)*xi) / 2  
    omega_complex = (xi.^2 + 1i .* (1 + 2*eta) .* xi) / 2;
    
    % Laplace exponent at complex argument
    lnL_main = laplace_exponent(omega_complex, p, alpha, T);
    
    % Final characteristic function 
    phi = exp(-1i .* xi .* lnL_eta + lnL_main);
end
