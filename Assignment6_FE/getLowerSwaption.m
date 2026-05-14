function [lowerValue, prices] = getLowerSwaption(dates, discounts, t0, callDates, paymentDates, K, a, sigma)
% Lower bound for Bermudan payer swaption: max European swaption.
% Each European priced in closed form via Jamshidian (1-factor Hull-White).
%
% NOTATION
%   B0(.) : initial discount factors (ZCB at t=0), via dates/discounts curve
%
% INPUT
%   dates, discounts - bootstrapped initial curve (provides B0)
%   t0               - settlement date (serial)
%   callDates        - vector of admissible exercise dates of the Bermudan
%                      (years 2..N-1 for non-call 2)
%   paymentDates     - ALL fixed-leg payment dates of the underlying swap
%                      (years 1..N)
%   K                - strike
%   a, sigma         - Hull-White parameters
%
% OUTPUT
%   lowerValue       - max of European swaption prices = lower bound
%   prices           - vector of European prices (one per call date)
%                      useful to know WHICH call date is optimal

    n = length(callDates);
    prices = zeros(n, 1);
    for j = 1:n
        T_alpha = callDates(j);
        T_pay   = paymentDates(paymentDates > T_alpha);
        if isempty(T_pay), continue; end
        prices(j) = europeanPayerSwaptionHW(dates, discounts, t0, ...
                                            T_alpha, T_pay, K, a, sigma);
    end
    lowerValue = max(prices);
end