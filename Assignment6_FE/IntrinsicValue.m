function V = IntrinsicValue(x, t, T, B0_t, B0_T, a, sigma, K)
% Payer swaption intrinsic value at exercise date t, per node.
%
% NOTATION
%   B0(.)    : initial discount factors (ZCB at t=0), observed today
%   D_t(.;x) : stochastic ZCB at tree node x, observed at exercise date t
%
% INPUTS
%   x       - (n_nodes x 1) OU state values at t
%   t       - scalar exercise date (year fraction from t0)
%   T       - row vector [T_{alpha+1}, ..., T_omega] of fixed-leg payment dates
%   B0_t    - B0(t) (scalar)
%   B0_T    - B0(T_k) row vector
%   a,sigma - Hull-White parameters
%   K       - swaption strike
%
% OUTPUT
%   V       - (n_nodes x 1) max(payer-swap value, 0) at each node

T = T(:)';

% Stochastic ZCB matrix D_t(T_k; x_j), one row per node, one column per payment
D_t = HJM_ZCB_vec(x, t, T, B0_t, B0_T, a, sigma);          % (n_nodes x n_pay)

% Fixed-leg year fractions: first accrual is t -> T_{alpha+1}
deltas = [T(1) - t, diff(T)];                              % (1 x n_pay)

% Annuity (BPV) at each node = sum_k delta_k * D_t(T_k; x)
BPV = D_t * deltas';                                       % (n_nodes x 1)

% Payer swap value at the node: floating leg minus fixed leg
V_swap = 1 - D_t(:, end) - K * BPV;

V = max(V_swap, 0);
end
