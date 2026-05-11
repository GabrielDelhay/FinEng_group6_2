function B = linearRateInterp(knownDates, knownDiscounts, settlementDate, t)
% LINEARRATEINTERP  Discount factor at one or many query dates by linear
%   interpolation of the zero rates.
%
%   B = linearRateInterp(knownDates, knownDiscounts, settlementDate, t)
%   takes the known curve (knownDates, knownDiscounts) and returns the
%   discount factor at every query date in t. Works transparently on a
%   scalar t or on an array of dates of any shape. The output B has the
%   same size and shape as t.
%
% INPUTS:
%   knownDates     : datenums of the known curve nodes
%   knownDiscounts : discount factors at those nodes
%   settlementDate : datenum of t0
%   t              : query datenum(s), scalar or array of any shape
%
% OUTPUT:
%   B              : discount factor(s), same shape as t

    inputSize = size(t);
    t_col     = t(:);                       % flatten for vectorized math

    % Drop the settlement node (B = 1 makes log(1)/0 indeterminate)
    valid  = knownDiscounts < 1 - 1e-14;
    kDates = knownDates(valid);
    kDisc  = knownDiscounts(valid);

    tau_known = (kDates - settlementDate) / 365;
    tau_t     = (t_col  - settlementDate) / 365;

    r_known = -log(kDisc) ./ tau_known;

    % Linear interpolation on the zero rates; flat extrapolation outside
    % the known knots.
    r_t = interp1(tau_known, r_known, tau_t, 'linear', 'extrap');

    B_col = exp(-r_t .* tau_t);

    % Force B = 1 exactly when the query coincides with the settlement,
    % so the result does not depend on interp1's extrapolation behaviour.
    B_col(tau_t == 0) = 1;

    B = reshape(B_col, inputSize);          % preserve the input shape
end % linearRateInterp