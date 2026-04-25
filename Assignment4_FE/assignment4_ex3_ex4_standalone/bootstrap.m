function [dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet)
% BOOTSTRAP  Bootstraps the Euribor 3M single-curve discount factor curve.
%
%   Uses three instrument types in sequence:
%     1. Deposits    – short end,  Act/360, simply compounded
%     2. STIR Futures– mid range,  fixed 0.25 year fraction (3M convention)
%     3. Swaps vs 3M – long end,   fixed leg annual 30/360 EU
%
% INPUTS:
%  datesSet: struct with settlementDate, deposDates, futuresDates, swapDates
%  ratesSet: struct with deposRates, futuresRates, swapRates
% OUTPUTS:
%   dates      – column vector of datenums (settlement first, then end dates)
%   discounts  – column vector of discount factors B(t0, T)
%   zeroRates  – column vector of zero rates (Act/365, continuously compounded)

settlementDate = datesSet.settlement;

% Mid-market rates (average of bid and ask)
midDepos   = mean(ratesSet.depos,   2);
midFutures = mean(ratesSet.futures, 2);
midSwaps   = mean(ratesSet.swaps,   2);

% Initialise curve at settlement: B(t0, t0) = 1
dates     = settlementDate;
discounts = 1.0;

%% DEPOSITS
%  B(t0, T) = 1 / (1 + r * delta)
%  Use first 3 deposits only (matches Python logic)

nDepos = length(datesSet.depos);
firstFutSettle = datesSet.futures(1,1);
for i = 1:nDepos
    T     = datesSet.depos(i); % end dates of depos
    delta = yearfrac(settlementDate, T, 2); %ACT/360
    B     = 1.0 / (1.0 + midDepos(i) * delta);
    [dates, discounts] = insertPoint(dates, discounts, T, B);
     if datesSet.depos(i) >= firstFutSettle
        break   % we stop when there is an overlap with futures (more liquid)
    end
end

%% STIR FUTURES
%  B(t0, T_i) = B(t0, T_{i-1}) * 1 / (1 + r_fwd * 0.25)
%  Chain on the last discount factor — no interpolation needed

nFut = min(7, size(datesSet.futures, 1));
for i = 1:nFut
    T1    = datesSet.futures(i, 1);   % period start (IMM date)
    T2   = datesSet.futures(i, 2);        % period end date
    B_T1 = linearRateInterp(dates, discounts, settlementDate, T1);
    delta = yearfrac(T1, T2, 2);          % Act/360

    B_T2 = B_T1 / (1.0 + midFutures(i) * delta);
    [dates, discounts] = insertPoint(dates, discounts, T2, B_T2);
end

% Ensure curve is sorted before swap bootstrap
[dates, sortIdx] = sort(dates);
discounts = discounts(sortIdx);

%% SWAPS vs Euribor 3M
% Single-curve bootstrap
% Fixed leg: annual payments, 30/360 European (basis 6)
% Bootstrap formula (n-th maturity):
% B(t0, T_n) = (1 - K * BPV_{1..n-1}) / (1 + K * delta_n)
% where  BPV = sum_{j=1}^{n-1}  delta_j * B(t0, T_j)
% B(t0, T_n) = (1 - K * BPV_{n-1}) / (1 + K * delta_n)
% BPV is maintained as a running accumulator across swap pillars

nSwaps   = length(datesSet.swaps);%number of swap at our disposal
% nFut is the number of fut we used
lastFutEnd = datesSet.futures(nFut, 2);

% Initialise BPV with the 1y swap (covered by futures, not bootstrapped)
firstSwapDate = datesSet.swaps(1);
yf_1 = yearfrac(settlementDate, firstSwapDate, 6);    % 30/360 EU
B_1  = linearRateInterp(dates, discounts, settlementDate, firstSwapDate);
BPV  = yf_1 * B_1;
prevDate = firstSwapDate;
for i = 1:nSwaps
    if datesSet.swaps(i) <= lastFutEnd
        continue   % skip pillars already covered by futures
    end

    T_n     = datesSet.swaps(i);
    K       = midSwaps(i);
    delta_n = yearfrac(prevDate, T_n, 6);              % 30/360 EU

    B_Tn = (1.0 - K * BPV) / (1.0 + K * delta_n);

    [dates, discounts] = insertPoint(dates, discounts, T_n, B_Tn);

    % Increment BPV for next iteration
    BPV      = BPV + delta_n * B_Tn;
    prevDate = T_n;
end

%% Zero rates
%  z(T) = -ln(B(t0,T)) / ((T - t0)/365)  Act/365 continuously compounded
T_ACT365  = (dates - settlementDate) / 365;
zeroRates = zeros(size(dates));

validIdx = T_ACT365 > 0;
zeroRates(validIdx) = -log(discounts(validIdx)) ./ T_ACT365(validIdx);

% Settlement point: zero rate undefined (T=0)
zeroRates(~validIdx) = NaN;

end % bootstrap


function [datesOut, discountsOut] = insertPoint(dates, discounts, t, B)
%INSERTPOINT  Insert or overwrite a point in the (sorted) curve.
idx = find(dates == t, 1);
if isempty(idx)
    datesOut     = [dates;     t];
    discountsOut = [discounts; B];
    [datesOut, sortIdx] = sort(datesOut);
    discountsOut = discountsOut(sortIdx);
else
    datesOut            = dates;
    discountsOut        = discounts;
    discountsOut(idx)   = B;
end
end % insertPoint