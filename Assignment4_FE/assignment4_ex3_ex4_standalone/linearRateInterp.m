function B = linearRateInterp(knownDates, knownDiscounts, settlementDate, t)
% Builds zero rates from the known curve (excluding settlement where B=1),
% linearly interpolates, flat-extrapolates, then converts back to DF.

% Exclude settlement (B = 1  =>  log(1)/0 = 0/0)
valid     = knownDiscounts < 1 - 1e-14;
kDates    = knownDates(valid);
kDisc     = knownDiscounts(valid);

tau_known = (kDates - settlementDate) / 365;
tau_t     = (t      - settlementDate) / 365;

r_known = -log(kDisc) ./ tau_known;

% Linear interpolation, flat extrapolation beyond last knot
r_t = interp1(tau_known, r_known, tau_t, 'linear', 'extrap');

B = exp(-r_t .* tau_t);
end % linearRateInterp