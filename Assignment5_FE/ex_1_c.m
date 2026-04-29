% Base setup and baseline pricing
formatDate = 'dd/mm/yyyy';
[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet);

BPV = 0.0001;
n_depos = 3;
n_futures = 7;
n_swaps = 50;
n_total = n_depos + n_futures + (n_swaps - 1);

% Pre-allocate price vector correctly (N x 1)
price = zeros(n_total, 1);

% Calculate baseline true price
true_price = ex2();

idx = 1; % Global index for the price vector

% 1. Bump Depos
for i = 1:n_depos
    ratesSet.depos(i) = ratesSet.depos(i) + BPV; % Bump
    [~, ~, ~] = bootstrap(datesSet, ratesSet);
    price(idx) = ex2();
    ratesSet.depos(i) = ratesSet.depos(i) - BPV; % Restore original value (CRITICAL)
    idx = idx + 1;
end

% 2. Bump Futures
for i = 1:n_futures
    ratesSet.futures(i) = ratesSet.futures(i) + BPV; 
    [~, ~, ~] = bootstrap(datesSet, ratesSet);
    price(idx) = ex2();
    ratesSet.futures(i) = ratesSet.futures(i) - BPV; 
    idx = idx + 1;
end

% 3. Bump Swaps (skip the first)
for i = 2:n_swaps
    ratesSet.swaps(i) = ratesSet.swaps(i) + BPV; % Fixed typo: was 'futures'
    [~, ~, ~] = bootstrap(datesSet, ratesSet);
    price(idx) = ex2();
    ratesSet.swaps(i) = ratesSet.swaps(i) - BPV; 
    idx = idx + 1;
end

% Vectorized delta calculation
delta = price - true_price;