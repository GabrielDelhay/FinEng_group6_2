% runAssignment6_Group6
%
%
clc; close all; clear all;

formatDate = 'dd/mm/yyyy';
[datesSet, ratesSet] = readExcelDataOS('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet);
t0 = dates(1);  % Settlement date: 19/02/2008


%% EXERCISE 1
% Calibration
load('eurostoxx_Poli.mat');
params_opt = calibration(cSelect, dates, zeroRates);   %sigma, kappa, eta from lab 4

% Contract Data
structData.S0 = 3200;                % Spot price at t0
structData.strike = 3200;            % Strike level for coupon calculation
structData.notional = 100e6;         % Principal Amount (100 MIO EUR)
structData.t0 = t0;                  % Settlement date: 19/02/2008
% Time structure
structData.T = [1, 2];               % Annual coupon observation dates (years)
structData.trigger = 0.06;           % Early Redemption Trigger Level (6%)
% Coupon logic
structData.couponYear1_Low  = 0.06;  % 6% if Stoxx50 < Strike
structData.couponYear1_High = 0;     % I'm actually not sure this is the right interpretation
structData.couponLastYear   = 0.02;  % Fixed 2% for the second year
% Floating Leg (What Bank XX pays to I.B.)
structData.spread = 0.0130;          % Euribor 3m + 1.30%
structData.freq   = 4;               % Quarterly payments (4 per year)

% Floating leg calculation (Euribor 3m + 1.30%)
% Schedule of Payment Dates
num_quarters = 8;
payment_dates = zeros(num_quarters, 1);
for k = 1:num_quarters
    % Increment by 3 months at each step from Start Date (t0)
    raw_date = addtodate(t0, 3*k, 'month');
    % Apply Modified Following convention using your provided function
    payment_dates(k) = adj_modfollow(raw_date);
end

% Calculate Year Fractions for Discouting and Accrual
t_discount = yearfrac(t0, payment_dates, 1); % Act/360 as per Annex 2 for discounting
% For the spread accrual: time between consecutive quarterly dates
period_starts = [t0; payment_dates(1:end-1)];
delta_pay = yearfrac(period_starts, payment_dates, 1); % Act/360 

% Filter out NaNs and interpolate zeroRates with spline
idxValid = ~isnan(zeroRates);     %first one is NaN: this solves it, maybe there is a better way
z_interp = interp1(dates(idxValid), zeroRates(idxValid), payment_dates, 'spline', 'extrap');
discounts_float = exp(-z_interp .* t_discount);

% Use your custom function to price the floating leg
% B_pay = discounts_float (interpolated with spline)
% delta_pay = tau_accrual (Act/360)
PV_floating_leg = price_floating_leg(structData.notional, discounts_float, delta_pay, structData.spread);
fprintf('Value of Floating Leg (via function): %.2f EUR\n', PV_floating_leg);


%% --- FIXED/CONDITIONAL LEG CALCULATION (NIG MODEL) ---

% 1. Get NIG parameters from Calibration
alpha = params_opt(1);
beta  = params_opt(2);
delta = params_opt(3);

% 2. Precise Dates and Discounts (Modified Following)
payment_dates = zeros(8, 1);
for k = 1:8
    raw_date = addtodate(t0, 3*k, 'month');
    payment_dates(k) = adj_modfollow(raw_date); % Using your function
end

% Annual coupon dates (T1 = end of 4th qtr, T2 = end of 8th qtr)
coupon_date1 = payment_dates(4); 
coupon_date2 = payment_dates(8);

% Interpolate discounts at exact dates
D1 = interp1(dates(idxValid), discounts(idxValid), coupon_date1, 'linear', 'extrap');
D2 = interp1(dates(idxValid), discounts(idxValid), coupon_date2, 'linear', 'extrap');

% Year fractions (Act/360 as per Annex 2)
T1_year = yearfrac(t0, coupon_date1, 2); 
T2_year = yearfrac(t0, coupon_date2, 2);

% 3. Pricing Digital Options via NIG
S0 = structData.S0;
K  = structData.strike;
r1_eff = -log(D1) / T1_year; 

drift_RN = delta * (sqrt(alpha^2 - beta^2) - sqrt(alpha^2 - (beta+1)^2));
phi_NIG = @(u, T, r) exp(1i*u*(log(S0) + (r - drift_RN)*T) + ...
    delta * T * (sqrt(alpha^2 - beta^2) - sqrt(alpha^2 - (beta - 1i*u).^2)));

integrand = @(u, T, r) real(exp(-1i*u*log(K)) .* phi_NIG(u, T, r) ./ (1i*u));
prob_S1_low = 0.5 - (1/pi) * integral(@(u) integrand(u, T1_year, r1_eff), 0, inf);

% 4. NPV and Early Redemption logic
% Year 1: 6% if Stoxx < Strike. Hits 6% trigger -> Early Redemption 
PV_coupon1 = structData.couponYear1_Low * structData.notional * prob_S1_low * D1;

% Survival: Swap continues only if Year 1 coupon was 0% (Stoxx >= Strike)
prob_survival = 1 - prob_S1_low; 

% Year 2: Fixed 2% paid only if swap survived T1 
PV_coupon2 = structData.couponLastYear * structData.notional * prob_survival * D2;

% Total Conditional Leg Value
PV_conditional_leg = PV_coupon1 + PV_coupon2;


fprintf('Value of Conditional Leg (NIG): %.2f EUR\n', PV_conditional_leg);

% NPV of swap: IB receives Conditional Leg, pays Floating Leg
NPV_Swap_IB = PV_conditional_leg - PV_floating_leg;
fprintf('Total NPV of the Swap for IB: %.2f EUR\n', NPV_Swap_IB);

% Compute upfront
Upfront_EUR = PV_conditional_leg - PV_floating_leg;
Upfront_Percentage = (Upfront_EUR / structData.notional) * 100;

fprintf('UPFRONT:%.2f EUR (%.4f%% of notional)\n', Upfront_EUR, Upfront_Percentage);