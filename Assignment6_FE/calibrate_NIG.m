function [p_opt, alpha, S0, divYield, strikes, smiles] = ...
         calibrate_NIG(cSelect, dates, zRates, fmt)
    strikes  = double(cSelect.strikes);
    S0       = double(cSelect.reference);
    divYield = double(cSelect.dividends);
    smiles   = double(cSelect.surface);
    T        = double(cSelect.maturity);
    alpha    = 0.5;
    matDate  = datenum('15/02/2008', fmt) + 365;
    r        = interp1(dates, zRates, matDate, 'linear', 'extrap');
    B        = exp(-r * T);
    F0       = S0 * exp((r - divYield) * T);
    C_mkt    = BS_call(F0, strikes, B, smiles, T);
    obj      = @(p) sum((nMV_call_FFT(p, alpha, F0, strikes, B, T) - C_mkt).^2);
    nonlcon  = @(p) deal(-p(3) - (1-alpha)/(p(2)*p(1)^2), []);
    opts     = optimoptions('fmincon','Display','off','MaxIterations',5000);
    p_opt    = fmincon(obj,[0.20,2.0,4.0],[],[],[],[],[1e-4,1e-4,-50],[2,10,50],nonlcon,opts);
end