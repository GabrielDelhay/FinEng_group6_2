function cal = getCalendar(N, t0, nStepsPerYear)
% Build the calendar for the Bermudan swaption pricer.
%
% INPUT
%   N             - underlying swap maturity in years
%   t0            - settlement serial date
%   nStepsPerYear - tree sub-steps per year
%
% OUTPUT
%   cal - struct with two fields:
%       cal.tree.dates    (M+1 x 1) serial dates of every tree node, t0 included
%       cal.tree.yf       (M+1 x 1) year fractions from t0 (ACT/365)
%       cal.tree.M        scalar number of intervals (= M)
%
%       cal.anchor.dates  (N x 1)   yearly swap anchors, modified following
%       cal.anchor.yf     (N x 1)   year fractions of anchors
%       cal.anchor.idx    (N x 1)   indices of anchors within cal.tree.dates
%                                   (closest-node snapping)

    M = N * nStepsPerYear;

    %% Tree grid (uniform in year fraction, includes t0)
    yfGrid    = (0:M)' / nStepsPerYear;                 % yfGrid(1) = 0
    nodeDates = t0 + round(yfGrid * 365.25);

    cal.tree.dates = nodeDates;
    cal.tree.yf    = yearfrac(t0, nodeDates, 3);
    cal.tree.M     = M;

    %% Annual swap anchors (modified following business day)
    anchorDates = arrayfun(@(i) addtodate(t0, i, 'year'), 1:N)';
    anchorDates = busdate(anchorDates - 1, 'modifiedfollow');

    %% Snap each anchor to the closest tree node
    anchorIdx = zeros(N, 1);
    for k = 1:N
        [~, anchorIdx(k)] = min(abs(nodeDates - anchorDates(k)));
    end

    cal.anchor.dates = anchorDates;
    cal.anchor.yf    = yearfrac(t0, anchorDates, 3);
    cal.anchor.idx   = anchorIdx;
end