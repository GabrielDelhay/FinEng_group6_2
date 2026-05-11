function [p_opt] = calibration(cSelect, dates, zeroRates)
% CALIBRATION - Calibrates the NIG model parameters (nMV parametrization)
%
% USAGE:
%   [p_opt] = calibration(cSelect, dates, zeroRates)
%
% INPUTS:
%   cSelect   - Struct containing market data:
%                 .strikes   : Vector of option strikes
%                 .maturity  : Time to maturity (years)
%                 .reference : Spot price (S0)
%                 .dividends : Dividend yield
%                 .surface   : Vector of market implied volatilities
%   dates     - Vector of dates for the bootstrap curve
%   zeroRates - Vector of zero-coupon rates corresponding to dates
%
% OUTPUTS:
%   p_opt     - 1x3 vector of optimal parameters [sigma, kappa, eta]
%               Note: These define the NIG model in its nMV representation.
%               The alpha (tail heaviness) is implicitly determined by these.

    %% 1. Data Extraction and Preparation
    strikes = double(cSelect.strikes); 
    T = double(cSelect.maturity);
    S0 = double(cSelect.reference); 
    divYield = double(cSelect.dividends);
    smiles_mkt = double(cSelect.surface);
    
    % Define date format for consistency
    formatDate = 'dd/mm/yyyy';
    
    % Interpolate the risk-free rate for the specific maturity
    % (Assumes 1-year maturity relative to the reference date)
    maturityDate = datenum('15/02/2008', formatDate) + 365;
    r = interp1(dates, zeroRates, maturityDate, 'linear', 'extrap');
    
    % Discount factor and Forward price calculation
    B = exp(-r * T);
    F0 = S0 * exp((r - divYield) * T);
    
    % Convert market Implied Volatilities into Call prices using Black-Scholes
    C_mkt = BS_call(F0, strikes, B, smiles_mkt, T);
    
    %% 2. Optimization Setup
    % Note: alpha_fft is the damping factor for the FFT price integration.
    % It is NOT the NIG alpha parameter.
    alpha_fft = 2/3; 
    
    % Initial guess for parameters: [sigma, kappa, eta]
    params0 = [0.20, 2.0, 4.0]; 
    
    % Constraints: sigma > 0, kappa > 0, eta is unbounded within limits
    lb = [1e-4, 1e-4, -50.0]; 
    ub = [2.00, 10.0, 50.0];
    
    % Nonlinear constraint to ensure the existence of the characteristic function
    nonlcon = @(p) deal(-p(3) - (1-alpha_fft)/(p(2)*p(1)^2), []);
    
    % Objective function: Sum of Squared Errors (SSE) between model and market prices
    obj = @(p) sum((nMV_call_FFT(p, alpha_fft, F0, strikes, B, T) - C_mkt).^2);
    
    % Solver options
    options = optimoptions('fmincon', 'Display', 'off', ...
        'MaxIterations', 5000, ...
        'OptimalityTolerance', 1e-10, ...
        'StepTolerance', 1e-10);
    
    %% 3. Run Calibration
    % fmincon finds the [sigma, kappa, eta] that minimize the objective function
    [p_opt, ~] = fmincon(obj, params0, [], [], [], [], lb, ub, nonlcon, options);

end