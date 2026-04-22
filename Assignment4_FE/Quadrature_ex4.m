function C_quad = Quadrature_ex4(x_vec, sigma, kappa, eta, T, alpha, B, F0)
C_quad = zeros(length(x_vec),1);
for j = 1:length(x_vec)
    x = x_vec(j);
    Integrand = @(u) real(lewisIntegrand_ex4(u, x, sigma, kappa, eta, T, alpha));
    I_q = 2 * quadgk(Integrand, 0, Inf); % Integral is equal to the double of the integral of the real part.
    C_quad(j) = B * F0 * (1 - exp(-x/2) / (2*pi) * I_q);
end