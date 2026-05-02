%% Exotic Cap - Calibration & Pricing under Bond Market Model (BMM)
clc; close all; clear all;

%% 1. Bootstrap
formatDate = 'dd/mm/yyyy';
[datesSet, ratesSet] = readExcelDataOS('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, ~] = bootstrap(datesSet, ratesSet);
t0 = datesSet.settlement;

%% 2. Market data
N          = 16;
resetDates = buildResetDates(t0, N);
delta      = yearfrac(resetDates(1:end-1), resetDates(2:end), 3);
B0_T       = linearRateInterp(dates, discounts, t0, resetDates);
L0         = (B0_T(1:end-1) ./ B0_T(2:end) - 1) ./ delta;

%% 3. Flat ATM vols from market table
maturitiesYears = [1, 2, 3, 4];
strikes_table   = [1.50,1.75,2.00,2.25,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,10.0] / 100;
vols_table = [
    0.140,0.130,0.129,0.121,0.133,0.138,0.144,0.150,0.172,0.191,0.202,0.216,0.239;
    0.224,0.197,0.175,0.180,0.192,0.204,0.210,0.214,0.223,0.236,0.249,0.261,0.281;
    0.238,0.217,0.200,0.198,0.203,0.205,0.208,0.214,0.229,0.243,0.256,0.267,0.282;
    0.242,0.224,0.209,0.204,0.204,0.202,0.202,0.205,0.217,0.229,0.240,0.250,0.266;
];
flatVolsATM_annual = getATMFlatVols(B0_T, delta, maturitiesYears, strikes_table, vols_table);

%% 4. BMM spot vol calibration
nu = spotvolbootstrap(maturitiesYears, flatVolsATM_annual, t0, resetDates, B0_T, delta, L0);

%% 5. Correlation matrix rho_{ij} = exp(-lambda * |T_i - T_j|)
lambda = 0.1;
tGrid  = yearfrac(t0, resetDates(2:end), 3);
RHO    = exp(-lambda * abs(tGrid' - tGrid));
C      = chol(RHO, 'lower');

%% 6. Monte Carlo pricing
Nsim    = 100000;
payoffs = priceMC_ExoticCap(nu, B0_T, delta, L0, C, RHO, N, Nsim);

price   = mean(payoffs);
std_err = std(payoffs) / sqrt(Nsim);
fprintf('Exotic Cap price (BMM): %.8f\n', price);
fprintf('Standard error:         %.8f\n', std_err);
fprintf('95%% CI: [%.8f, %.8f]\n', price - 1.96*std_err, price + 1.96*std_err);
