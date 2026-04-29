%% run1b -- Pricing of the structured bond 

clear; clc;

% Run point 1.a to get the calibrated LMM spot-vol surface and all the
% supporting quantities (cap_dates, fwd_rates, delta_fwd, tau_expiry,
% all_B, spot_vols, strikes, cap_maturity_idx, t0, ...).
EX_1;

%% Bond parameters from the termsheet

N        = 50e6;       % notional, 50 MIO EUR
spread = 0.0200;       % Party A pays Euribor 3m + 2.00%

bond.first_cpn = 0.04;
bond.spread1 = 0.0100;   bond.K1 = 0.0420;   bond.c_dig1 = 0.0070;  % Regime 1
bond.spread2 = 0.0120;   bond.K2 = 0.0470;   bond.c_dig2 = 0.0100;  % Regime 2 
bond.spread3 = 0.0130;   bond.K3 = 0.0540;   bond.c_dig3 = 0.0110;  % Regime 3

%% Bond schedule
%
% Convention used everywhere below: period i pays at T_i (= cap_dates(i))
% and resets at T_{i-1}. Period 1 is the first fixed coupon (4%) paid at
% Start + 3m; periods 2..40 are the structured coupons.
%
% In EX_1 cap_maturity_idx was decremented by 1 for the bootstrap loop,
% so we add 1 to recover the index of the period whose payment date
% matches the regime boundary (3y, 6y, 10y).

i_3y  = cap_maturity_idx(3)  + 1;     % period paying at 3y  (i = 12)
i_6y  = cap_maturity_idx(6)  + 1;     % period paying at 6y  (i = 24)
i_10y = cap_maturity_idx(10) + 1;     % period paying at 10y (i = 40)

schedule.B_pay     = all_B(2 : i_10y+1);     % B(t0, T_i) for i=1..i_10y
schedule.delta_pay = delta_fwd(1 : i_10y);   % delta(T_{i-1}, T_i)
schedule.tau_reset = tau_expiry(1 : i_10y);  % (T_{i-1}-t0)/365
schedule.fwd_rates = fwd_rates(1 : i_10y);   % L_{i-1}(t0)
schedule.i_first   = 2;                      % first STRUCTURED coupon
schedule.i_3y      = i_3y;
schedule.i_6y      = i_6y;
schedule.i_10y     = i_10y;

%% Market data (vols + strikes)

market.spot_vols   = spot_vols(1 : i_10y, :);   % aligned with schedule
market.strikes_mkt = strikes / 100;             % decimals (in EX_1 strikes are %)

%% Pricing

NPV_A = price_floating_leg(N, schedule.B_pay, schedule.delta_pay, spread);
NPV_B = price_coupon_leg(N, schedule, market, bond);

X = (NPV_A - NPV_B) / N;

%% Output

fprintf('\n=== Pricing of the structured bond (point 1.b) ===\n');
fprintf('Regime 3 variant       : %s\n', bond.regime3_type);
fprintf('Notional               : %.0f EUR\n', N);
fprintf('NPV Party A (floating) : %14.2f EUR\n', NPV_A);
fprintf('NPV Party B (coupons)  : %14.2f EUR\n', NPV_B);
fprintf('Upfront X*N            : %14.2f EUR\n', X*N);
fprintf('Upfront X              : %10.4f %%\n', X*100);