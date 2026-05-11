% runAssignment6_Group6
%
%
clc; close all; clear all;


formatDate = 'dd/mm/yyyy';
[datesSet, ratesSet] = readExcelDataOS('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet);
t0 = dates(1);  % Settlement date: 19/02/2008


%% EXERCISE 1
