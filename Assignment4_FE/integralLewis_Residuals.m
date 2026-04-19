function C = integralLewis_Residuals(x, p_plus, p_minus, mu)

i = 1i;

if x > 0
    % pôle dans le demi-plan supérieur
    u = 1i*(p_minus - 0.5);

elseif x < 0
    % pôle dans le demi-plan inférieur
    u = -1i*(p_plus + 0.5);

else
    % cas limite
    Cp = integralLewis_Residuals(1e-8, p_plus, p_minus, mu);
    Cm = integralLewis_Residuals(-1e-8, p_plus, p_minus, mu);
    C = 0.5*(Cp + Cm);
    return;
end

v = -u - 1i/2;

phi = exp(1i * mu * v) ./ ...
    ((1 - 1i * v / p_plus) .* (1 + 1i * v / p_minus));

% dérivée correcte du facteur nul uniquement
if x > 0
    d = (1 - 1i*v/p_plus) * (1i/p_minus);
else
    d = (1 + 1i*v/p_minus) * (-1i/p_plus);
end

Res = phi * exp(-1i*x*u) / d;

C = 2*pi*1i * Res;

if x < 0
    C = -C;
end

C = real(C);

end

