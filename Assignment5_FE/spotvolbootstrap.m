function nu = spotvolbootstrap(maturitiesYears, flatVolsATM_annual, t0, ...
                               resetDates, B0_T, delta, L0)
% Bootstraps BMM spot vols nu(1..15) from flat Black cap vols by bucket
% (1Y,2Y,3Y,4Y), matching market cap prices. Within each annual bucket,
% nu is linearly interpolated between the boundary vols.

% Inputs: maturitiesYears, 
%         flatVolsATM_annual, 
%         t0, 
%         resetDates,
%         B0_T, 
%         delta, 
%         L0

    % Number of forward rates / caplets
    nCaplets = 15;

    % Output vector
    nu = zeros(1, nCaplets);

    % Market cap prices from ATM Black formula
    %
    % For ATM caplets:
    %   K_i = L0(i+1)
    % therefore:
    %   d1 = +0.5 * sigma * sqrt(T)
    %   d2 = -0.5 * sigma * sqrt(T)
    %
    % Cap maturities: 1Y, 2Y, 3Y, 4Y


    capPriceMarket = zeros(1,4);

    for m = 1:4

        nQ = 4 * m;

        % Caplet indices contributing to the cap
        idx = 1:(nQ-1);

        % Fixing times T_{i+1}
        Ti = yearfrac(t0, resetDates(idx+2), 3);

        % Flat ATM volatility associated with the cap maturity
        capMat = yearfrac(t0, resetDates(nQ+1), 3);

        fv = interp1(maturitiesYears, ...
                     flatVolsATM_annual, ...
                     capMat, ...
                     'linear', ...
                     'extrap');

        % ATM Black quantities
        d1 = 0.5 * fv .* sqrt(Ti);
        d2 = -d1;

        % Vectorized cap price
        capPriceMarket(m) = sum( ...
            B0_T(idx+2) .* ...
            delta(idx+1) .* ...
            L0(idx+1) .* ...
            (normcdf(d1) - normcdf(d2)) ...
        );

    end

    % Bootstrap first bucket (1Y)
    %
    % Constant instantaneous volatility over the first bucket:
    %   nu(1:3) = constant

    target_1Y = capPriceMarket(1);

    idx1 = 1:3;

    f1 = @(nu_cst) ...
        sum( ...
            arrayfun(@(i) ...
                capletBMM_strike_ATM( ...
                    nu_cst, ...
                    B0_T(i+2), ...
                    delta(i+1), ...
                    L0(i+1), ...
                    yearfrac(t0, resetDates(i+2), 3)), ...
                idx1)) ...
        - target_1Y;

    nu_1Y = fzero(f1, 0.002);

    nu(idx1) = nu_1Y;

    sigma_alpha = nu_1Y;

    % Price already explained by calibrated buckets
    cap_price_already = sum( ...
        arrayfun(@(i) ...
            capletBMM_strike_ATM( ...
                nu(i), ...
                B0_T(i+2), ...
                delta(i+1), ...
                L0(i+1), ...
                yearfrac(t0, resetDates(i+2), 3)), ...
            idx1));

    % Bootstrap remaining buckets: 2Y, 3Y, 4Y
    %
    % Volatility is linearly interpolated between:
    %   sigma_alpha -> sigma_beta

    for k = 2:4

        % Bucket boundaries
        i_alpha = (k-1) * 4;
        i_beta  = min(k * 4 - 1, nCaplets);

        idx = i_alpha:i_beta;

        T_alpha = resetDates(i_alpha + 1);
        T_beta  = resetDates(i_beta  + 1);

        % Objective function:
        % previously calibrated price
        % + contribution of current bucket
        f = @(sigma_beta) ...
            cap_price_already + ...
            priceBMM_segment( ...
                sigma_beta, ...
                sigma_alpha, ...
                i_alpha, ...
                i_beta, ...
                B0_T, ...
                delta, ...
                L0, ...
                resetDates, ...
                t0, ...
                T_alpha, ...
                T_beta) ...
            - capPriceMarket(k);

        % Initial guess:
        % approximate conversion LMM -> BMM
        weights = (delta(idx) .* L0(idx+1)) ./ ...
                  (1 + delta(idx) .* L0(idx+1));

        nu_init = flatVolsATM_annual(k) * mean(weights);

        sigma_beta = fzero(f, nu_init);

        % Linear interpolation inside the bucket
        w = (resetDates(idx+1) - T_alpha) ./ (T_beta - T_alpha);

        nu(idx) = sigma_alpha + w .* (sigma_beta - sigma_alpha);

        % Prepare next bucket
        sigma_alpha       = sigma_beta;
        cap_price_already = capPriceMarket(k);

    end

    %% Display calibrated spot volatilities

    fprintf('\nBootstrapped instantaneous volatilities (BMM):\n');
    disp(nu);

end
