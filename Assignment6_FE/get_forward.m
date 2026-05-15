function [F, tau, x_K] = get_forward(t0, T_reset, S0, divYield, dates, zeroRates, K)
% Computes the forward price, time to reset, and log-moneyness at a reset date.
%
% INPUTS:
    % t0        : pricing date (MATLAB serial date)
    % T_reset   : reset date (MATLAB serial date)
    % S0        : spot price
    % divYield  : continuous dividend yield
    % dates     : bootstrapped dates vector
    % zeroRates : continuously compounded zero rates matching dates
    % K         : strike price
%
% OUTPUTS:
    % F         : forward price F(t0, T_reset) = S0 * exp((r - q) * tau)
    % tau       : yearfrac(t0, T_reset) Act/365
    % x_K       : log-moneyness log(K/F)
%
    tau = yearfrac(t0, T_reset, 3);
    r   = interp1(dates, zeroRates, T_reset, 'linear', 'extrap');
    F   = S0 * exp((r - divYield) * tau);
    x_K = log(K / F);
end