function [p_inner, p_bot, p_top] = treeProbabilities(a, dt, lmax)
% Trinomial tree transition probabilities (Hull-White, Baviera convention).
% Branching patterns:
%   - inner  (|l| < lmax) : children at (l+1, l, l-1)  -- pattern A
%   - bottom (l = -lmax)  : children at (l+2, l+1, l)  -- pattern B
%   - top    (l = +lmax)  : children at (l, l-1, l-2)  -- pattern C
% In every output, columns are [p_to_up_target, p_to_mid_target, p_to_down_target].
%
% Inputs:
%   a    - mean reversion speed
%   dt   - tree time step
%   lmax - truncation level (positive integer)
%
% Outputs:
%   p_inner - (2*lmax-1) x 3 matrix; row r corresponds to level l = r - lmax
%   p_bot   - 1 x 3 row vector for l = -lmax
%   p_top   - 1 x 3 row vector for l = +lmax

mu_hat = 1 - exp(-a*dt);    % exact OU mean-revert factor per step

% ---- Pattern A: inner levels ----
l = (lmax-1 : -1 : -lmax+1).';
y = l * mu_hat;
p_inner = [ 0.5*(1/3 - y + y.^2), ...
    2/3 - y.^2,           ...
    0.5*(1/3 + y + y.^2) ];

% ---- Pattern B: bottom (l = -lmax) ----
yb = -lmax * mu_hat;
p_bot = [ 0.5*(1/3 + yb + yb^2),  ...
    -1/3 - 2*yb - yb^2,     ...
    0.5*(7/3 + 3*yb + yb^2) ];

% ---- Pattern C: top (l = +lmax) ----
yt = lmax * mu_hat;
p_top = [ 0.5*(7/3 - 3*yt + yt^2), ...
    -1/3 + 2*yt - yt^2,      ...
    0.5*(1/3 - yt + yt^2)   ];
end