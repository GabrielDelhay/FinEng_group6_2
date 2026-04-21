%  Compare the price of a cash-or-nothing digital call option computed:
%    (1) Under the Black model (flat volatility = ATM vol, no smile)
%    (2) Taking into account the implied volatility smile (call-spread replication)
%
%  CONTRACT SPECS:
%    - Strike       : ATM forward
%    - Expiry       : 1 year (Act/365)

clear; clc; close all;

% Load the dataset (struct 'cSelect' contains all market info)
load('eurostoxx_Poli.mat');

strikes = double(cSelect.strikes);   % [31x1] vector of strikes
T = double(cSelect.maturity);   % Time to maturity in years 
S0 = double(cSelect.reference);   % ATM Spot price = 3771.06 
divYield = double(cSelect.dividends);   % Continuous dividend yield
volSmile = double(cSelect.surface);   % [31x1] implied volatility smile

formatDate = 'dd/mm/yyyy';
valuationdate = '15/02/2008';
valuationDate = datenum(valuationdate, formatDate);
maturityDate = valuationDate + 365;   % Act/365, T=1y
[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet);
r = interp1(dates, zeroRates, maturityDate, "linear", "extrap");

B = exp(-r * T);         % Discount factor B(0, T) = exp(-r*T)

%  ATM Forward price
F0 = S0 * exp((r - divYield) * T);

%Contract parameters
Notional      = 10e6;
digitalPayoff = 0.05 * Notional;          % 5% of Notional (500,000 EUR) (paid if S_T > K)
K_digital = F0;      % strike of the digital

% Interpolate smile at the ATM strike K = F0 to obtain ATM implied volatility
sigma_ATM = interp1(strikes, volSmile, K_digital, 'spline');  %we choose spline because we need differentiability

% The Black-Scholes price of a cash-or-nothing digital CALL is: (NO SMILE)
% Digital_Black = B * cashAmount * N(d2)

d1_black = ( log(F0 / K_digital) + 0.5 * sigma_ATM^2 * T ) / (sigma_ATM * sqrt(T));
d2_black = d1_black - sigma_ATM * sqrt(T);

% Price of digital = discounted expected payoff under risk-neutral measure
price_Black = B * digitalPayoff * normcdf(d2_black)

%  SMILE-ADJUSTED PRICING (call-spread replication)
%
%  A cash-or-nothing digital can be replicated as the limit of a call spread:
%
%    Digital(K) = lim_{eps -> 0}  [ C(K - eps) - C(K + eps) ] / (2*eps)
%
%  In discrete form with finite epsilon:
%    Digital(K) ≈ [ C(K - eps) - C(K + eps) ] / (2*eps)
%
%  Each call price C(K +/- eps) is computed using the Black formula with
%  the SMILE-interpolated implied volatility at that specific strike.
%  This naturally incorporates the slope of the smile (the "vanna" effect).
%
%  We choose eps small enough to be accurate but not so small as to cause
%  numerical instability in the interpolation.
% =========================================================================

eps = 1.0;    % 1 index point (small relative to the ~3800 level of the index)

% Strikes slightly below and above the digital strike
K_lo = K_digital - eps;
K_hi = K_digital + eps;

% Interpolate smile at the shifted strikes (cubic spline)
sigma_lo = interp1(strikes, volSmile, K_lo, 'spline');
sigma_hi = interp1(strikes, volSmile, K_hi, 'spline');

% Compute Black call prices at the two shifted strikes using smile vols
C_lo = blackCallPrice(F0, K_lo, sigma_lo, T, B);
C_hi = blackCallPrice(F0, K_hi, sigma_hi, T, B);

% Approximate the digital price via the call spread
% Note: we multiply by cashAmount to get the actual EUR payoff
price_Smile = digitalPayoff * (C_lo - C_hi) / (2 * eps)


%% =========================================================================
%  SECTION 8: COMPARISON AND DISCUSSION
% =========================================================================

diff_EUR  = price_Smile - price_Black;
diff_perc = 100 * diff_EUR / price_Black;


%% =========================================================================
%  SECTION 9: PLOT THE IMPLIED VOLATILITY SMILE
% =========================================================================

% Extended fine grid for smooth smile plot
K_fine = linspace(strikes(1), strikes(end), 500);
vol_fine = interp1(strikes, volSmile, K_fine, 'spline');

figure('Name', 'Implied Volatility Smile - Eurostoxx 50 (15 Feb 2008)');
plot(strikes, volSmile*100, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6); hold on;
plot(K_fine, vol_fine*100, 'b-', 'LineWidth', 1.5);
xline(F0, 'r--', 'LineWidth', 1.5, 'Label', sprintf('ATM F_0 = %.1f', F0));
xline(K_digital, 'k:', 'LineWidth', 1.2);
ylabel('Implied Volatility (%)');
xlabel('Strike K');
title('Eurostoxx 50 Implied Volatility Smile - 1Y, 15 Feb 2008');
legend('Market vols', 'Spline interpolation', 'ATM Forward', 'Location', 'NorthEast');
grid on;


%% =========================================================================
%  LOCAL FUNCTION: Black Call Price (normalized, i.e. for a unit payoff)
%
%  Computes the Black-Scholes call price given:
%    F0    : ATM forward price
%    K     : strike
%    sigma : implied volatility at this strike
%    T     : time to maturity
%    B     : discount factor B(0,T)
%
%  Returns the call price per unit of notional (i.e. as a fraction).
%  To get the EUR price: multiply by Notional.
% =========================================================================

function C = blackCallPrice(F0, K, sigma, T, B)
    % d1 and d2 for Black formula (forward-based, no drift term needed)
    d1 = ( log(F0 / K) + 0.5 * sigma^2 * T ) / (sigma * sqrt(T));
    d2 = d1 - sigma * sqrt(T);
    
    % Black call price (discounted, forward-based)
    C = B * ( F0 * normcdf(d1) - K * normcdf(d2) );
    
    % Normalize by K so the call spread gives a probability-like quantity
    % Actually: keep in index points, the caller divides by (2*eps) and
    % multiplies by cashAmount -- so C is in index points here.
end