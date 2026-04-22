function [F0, B, r, T] = get_datas()
%% Bootstrap
formatDate = 'dd/mm/yyyy';
maturity = datenum('15-Feb-2009');
[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet);

%% Exercise 4
%% data
load('eurostoxx_Poli.mat');

S0       = cSelect.reference;          % spot
T        = cSelect.maturity;           % 1 year
q        = cSelect.dividends;          % continuous dividend yield

% Discount factor from bootstrap
r = interp1(dates, zeroRates, maturity);                      
B = exp(-r * T);
F0 = S0 * exp((r - q) * T);      % ATM forward
end