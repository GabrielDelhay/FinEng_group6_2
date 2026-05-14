clc; close all; clear all;

formatDate = 'dd/mm/yyyy';
[datesSet, ratesSet] = readExcelDataOS('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, ~] = bootstrap(datesSet, ratesSet);
t0 = dates(1);

% Parametri modello
sigma = 0.8e-2;
a     = 0.11;
K     = 0.05;
N     = 10;

% Tree settings
nStepsPerYear = 365;
dt    = 1/nStepsPerYear;
dx    = sigma * sqrt(3*dt);
muHat = 1 - exp(-a*dt);
lMax  = ceil((1 - sqrt(2/3))/muHat);
x     = dx * (lMax:-1:-lMax)';

%% Calendar (struct)
cal = getCalendar(N, t0, nStepsPerYear);
M   = cal.tree.M;

% Pre-compute initial DFs at each tree node
B0_node = arrayfun(@(d) linearRateInterp(dates, discounts, t0, d), cal.tree.dates);

% Pre-compute deltaX matrices and tree probabilities
[deltaX_top, deltaX_inner, deltaX_bot] = builddeltaX(dx, lMax);
[p_inner, p_bot, p_top] = treeProbabilities(a, dt, lMax);

%% Anchor subset for the swap (drop the non-call year 1, keep maturity)
exerciseYF   = cal.anchor.yf(2:end);                 % years 2..10
exercise_idx = cal.anchor.idx(2:end);

%% 2.A — Bermudan price via tree
V_terminal = zeros(2*lMax+1, 1);

priceSwaption = bermudanRollback(V_terminal, x, ...
    p_inner, p_bot, p_top, ...
    dt, dx, ...
    cal.tree.yf, B0_node, ...
    exerciseYF, exercise_idx, ...
    a, sigma, K, lMax, 1);

fprintf('Bermudan payer swaption price : %.6f\n', priceSwaption);

%% 2.B — Tree sanity check: reprice B(0, T_N)
V_terminal = ones(2*lMax+1, 1);

B0_10_tree = bermudanRollback(V_terminal, x, ...
    p_inner, p_bot, p_top, ...
    dt, dx, ...
    cal.tree.yf, B0_node, ...
    exerciseYF, exercise_idx, ...
    a, sigma, K, lMax, 0);                           % allowExercise = 0!

B0_10 = linearRateInterp(dates, discounts, t0, cal.tree.dates(end));
fprintf('B(0,T_N) tree vs curve        : %.2e\n', B0_10 - B0_10_tree);

%% 2.C — Bounds
exec = cal.anchor.dates(2:end);                      % T_2..T_10

upperBound = getUpperSwaption(dates, discounts, t0, exec, K, a, sigma);
fprintf('Upper bound (cap)             : %.6f\n', upperBound);

callDates    = cal.anchor.dates(2:end-1);            % T_2..T_9
paymentDates = cal.anchor.dates;                     % T_1..T_10

lowerBound = getLowerSwaption(dates, discounts, t0, ...
    callDates, paymentDates, K, a, sigma);
fprintf('Lower bound (max Eur Jamsh.)  : %.6f\n', lowerBound);

%% Recap
fprintf('\n=== Summary ===\n');
fprintf('Lower    : %.6f\n', lowerBound);
fprintf('Bermudan : %.6f\n', priceSwaption);
fprintf('Upper    : %.6f\n', upperBound);
fprintf('In band? : %s\n', string(priceSwaption >= lowerBound && priceSwaption <= upperBound));