function V0 = bermudanRollback(V_terminal, x, ...
    p_inner, p_bot, p_top, ...
    dt, dx, ...
    yfGrid, B0_node, ...
    exerciseYF, exercise_idx, ...
    a, sigma, K, lMax, allowExercise)
% Backward induction on the HW trinomial tree for a Bermudan payer swaption.
% Storage convention: x in DESCENDING order (x(1) = +lmax*deltaX, x(end) = -lmax*deltaX).
%
% NOTATION
%   B0(.)    : initial discount factors (ZCB at t=0) from the market curve
%   D_curr   : stochastic one-period ZCB at the current node (from HJM_ZCB_vec)
%   D_inner/top/bot : realized stochastic discount per (node, branch)
%
% INPUTS
%   V_terminal    - terminal payoff vector at yfGrid(end), nNodes x 1
%   x             - OU node values, descending, nNodes x 1
%   p_inner       - (nNodes-2) x 3 inner transition probabilities (pattern A)
%   p_bot         - 1 x 3 boundary probabilities at x(end), pattern B
%   p_top         - 1 x 3 boundary probabilities at x(1),   pattern C
%   dt, dx        - tree step in time and in x
%   yfGrid        - 1 x M year fractions from t0 (yfGrid(1) = 0)
%   B0_node       - 1 x M initial discount factors B0(.) at each grid date
%   exerciseYF    - 1 x N year fractions of yearly anchors [T_1, ..., T_N]
%   exercise_idx  - 1 x N indices of those anchors within yfGrid
%   a, sigma      - Hull-White parameters
%   K             - payer swaption strike
%   allowExercise - (optional, default true) toggle for early exercise
%
% OUTPUT
%   V0            - scalar option value at t = 0


V_old = V_terminal(:);
M     = numel(yfGrid);

% Per-branch signed displacements in x (built once)
[dX_top, dX_inner, dX_bot] = builddeltaX(dx, lMax);



for idx = M-1 : -1 : 1
    t_curr = yfGrid(idx);
    t_next = yfGrid(idx+1);

    % 1) one-period stochastic ZCB at each node: D_curr(j) = D_{t_curr}(t_next; x_j)
    D_curr = HJM_ZCB_vec(x, t_curr, t_next, B0_node(idx), B0_node(idx+1), a, sigma);

    % 2) realized stochastic discounts per (node, branch)
    D_inner = stochasticDiscount_fwd(x(2:end-1), dX_inner, a, sigma, dt, D_curr(2:end-1));
    D_top   = stochasticDiscount_fwd(x(1),       dX_top,   a, sigma, dt, D_curr(1));
    D_bot   = stochasticDiscount_fwd(x(end),     dX_bot,   a, sigma, dt, D_curr(end));

    % 3) continuation value
    V_new = zeros(size(V_old));
    V_new(2:end-1) = p_inner(:,1).*D_inner(:,1).*V_old(1:end-2) + ...
        p_inner(:,2).*D_inner(:,2).*V_old(2:end-1) + ...
        p_inner(:,3).*D_inner(:,3).*V_old(3:end);
    V_new(1)   = p_top(1)*D_top(1)*V_old(1)     + ...
        p_top(2)*D_top(2)*V_old(2)     + ...
        p_top(3)*D_top(3)*V_old(3);
    V_new(end) = p_bot(1)*D_bot(1)*V_old(end-2) + ...
        p_bot(2)*D_bot(2)*V_old(end-1) + ...
        p_bot(3)*D_bot(3)*V_old(end);
    V_old = V_new;

    % 4) early exercise (skip non-call date T_1 and maturity T_N)
    if allowExercise && ismember(idx, exercise_idx)
        k        = find(exercise_idx == idx);
        T_pay_yf = exerciseYF(k+1:end);
        B0_pay   = B0_node(exercise_idx(k+1:end));
        IV = IntrinsicValue(x, t_curr, T_pay_yf(:).', ...
            B0_node(idx), B0_pay(:).', a, sigma, K);
        V_old = max(V_old, IV);
    end
end

% At t=0 the tree collapses to the center node (x = 0): middle index
V0 = V_old((numel(V_old)+1)/2);
end