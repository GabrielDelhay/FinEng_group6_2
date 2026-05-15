function [T, T_reset, B, BPV] = get_schedule(t0, months, dates, discounts, adj)
% Builds payment date, reset date, discount factor and floating leg BPV
% for a single year-end of the certificate schedule.
% BPV covers the last 4 quarterly periods ending at T (Act/360).
%
% INPUTS:
    % t0        : pricing date (MATLAB serial date)
    % months    : maturity in months (12, 24, or 36)
    % dates     : bootstrapped dates vector
    % discounts : discount factors matching dates
    % adj       : date adjustment function handle (following convention)
%
% OUTPUTS:
    % T         : adjusted payment date at t0 + months
    % T_reset   : reset date = 2 business days before T
    % B         : discount factor B(t0, T)
    % BPV       : sum of delta_i * B(t0, q_i) over the 4 quarters of the year
%    
    T       = adj(datemnth(t0, months));
    T_reset = busdate(busdate(T,-1),-1);
    B       = linearRateInterp(dates, discounts, t0, T);
    q       = adj_modfollow(datemnth(t0, (3:3:months)'));
    B_q     = linearRateInterp(dates, discounts, t0, q);
    prev_q  = [t0; q(1:end-1)];
    delta_q = yearfrac(prev_q, q, 2);
    BPV     = sum(delta_q(end-3:end) .* B_q(end-3:end));
end

