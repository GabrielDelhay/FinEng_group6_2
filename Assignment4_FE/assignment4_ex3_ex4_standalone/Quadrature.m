function C_quad = Quadrature(x_vec, phi, B, F0)
% Quadrature - Call prices via Lewis formula evaluated with quadgk.
%
% Inputs:
%   x_vec - moneyness vector log(F0/K)
%   phi   - characteristic function handle (martingale measure)
%   B, F0 - discount factor and ATM forward
%
% Output:
%   C_quad - column vector of call prices, same length as x_vec

x_vec  = x_vec(:);
C_quad = zeros(length(x_vec), 1);

for j = 1:length(x_vec)
    x = x_vec(j);
    Integrand = @(u) real(lewisIntegrand(u, x, phi));
    % integrand is even under u -> -u (Hermitian symmetry on the real axis)
    I_q = 2 * quadgk(Integrand, 0, Inf);
    C_quad(j) = B * F0 * (1 - exp(-x/2) / (2*pi) * I_q);
end

end
