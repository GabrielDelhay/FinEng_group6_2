function upperValue = getUpperSwaption(dates, discounts, t0, executionDates, K, a, sigma)
% Upper bound for the Bermudan payer swaption via a strip of caplets (a cap).
%
% Rationale:
%   max(0, sum_k delta_k (L_k - K)) <= sum_k delta_k max(0, L_k - K)
% so any payer swaption is bounded above by the cap on the same period;
% the Bermudan (one-shot exercise) is in turn bounded by such a cap.
%
% NOTATION
%   B0(.) : initial discount factors (ZCB at t=0), via dates/discounts curve
%
% INPUT
%   dates, discounts - bootstrapped initial curve (provides B0)
%   t0               - settlement date (serial)
%   executionDates   - vector of yearly anchors defining the strip of caplets.
%                      Pair (executionDates(i), executionDates(i+1)) defines
%                      the i-th caplet: reset at i, payment at i+1.
%                      For non-call 2 on a 10y swap pass [T_2, T_3, ..., T_10]
%                      (9 dates -> 8 caplets).
%   K                - strike
%   a, sigma         - Hull-White parameters
%
% OUTPUT
%   upperValue       - cap price = upper bound for the Bermudan

    upperValue = 0;
    for i = 1:length(executionDates)-1
        T_reset = executionDates(i);
        T_pay   = executionDates(i+1);
        upperValue = upperValue + ...
            capletHW(dates, discounts, t0, T_reset, T_pay, K, a, sigma);
    end
end