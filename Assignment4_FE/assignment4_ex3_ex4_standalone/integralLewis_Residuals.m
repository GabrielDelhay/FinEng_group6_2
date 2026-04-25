function C = integralLewis_Residuals(x, p_plus, p_minus, mu)

i = 1i;

if x > -mu
    % x>-mu: the contour must be closed in the lower half-plane -> I = -2*pi*i * sum(Res_LHP)
    % Res 1: pole of u^2+0.25 in u = -i/2
    u1 = -i/2;
    v1 = -u1 - i/2;           % v1 = 0  =>  phi(0) = 1
    phi1 = exp(i*mu*v1) / ((1 - i*v1/p_plus) * (1 + i*v1/p_minus));
    Res1 = phi1 * exp(-i*x*u1) / (2*u1);   % = i*exp(-x/2)

    % Res 2: pole of phi in u = -i*(p_minus+0.5), v2 = i*p_minus
    u2 = -i*(p_minus + 0.5);
    v2 = -u2 - i/2;
    phi_num       = exp(i*mu*v2);
    non_vanishing = 1 - i*v2/p_plus;
    deriv         = -i/p_minus;            % d/du of (1 + i*v/p_minus)
    Res2 = phi_num * exp(-i*x*u2) / (non_vanishing * deriv * (u2^2 + 0.25));

    I = -2*pi*i * (Res1 + Res2);

elseif x < -mu
    % x<-mu: the contour must be closed in the upper half-plane -> I = +2*pi*i * sum(Res_UHP)
    % Res 1: pole of u^2+0.25 in u = +i/2
    u1 = i/2;
    v1 = -u1 - i/2;           % v1 = -i  =>  phi(-i) = 1 (condition martingale)
    phi1 = exp(i*mu*v1) / ((1 - i*v1/p_plus) * (1 + i*v1/p_minus));
    Res1 = phi1 * exp(-i*x*u1) / (2*u1);   % = -i*exp(x/2)

    % Res 2: pole of phi in u = +i*(p_plus-0.5), v2 = -i*p_plus
    u2 = i*(p_plus - 0.5);
    v2 = -u2 - i/2;
    phi_num       = exp(i*mu*v2);
    non_vanishing = 1 + i*v2/p_minus;
    deriv         = i/p_plus;              % d/du of (1 - i*v/p_plus)
    Res2 = phi_num * exp(-i*x*u2) / (non_vanishing * deriv * (u2^2 + 0.25));

    I = 2*pi*i * (Res1 + Res2);

else
    Cp = integralLewis_Residuals(1e-8, p_plus, p_minus, mu);
    Cm = integralLewis_Residuals(-1e-8, p_plus, p_minus, mu);
    I = 0.5*(Cp + Cm);
end

C = real(I);

end


