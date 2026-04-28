%% Exotic Cap - Calibration & Pricing under Bond Market Model (BMM)

%% 1. Bootstrap
formatDate = 'dd/mm/yyyy';
[datesSet, ratesSet] = readExcelDataOS('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet);
t0 = datesSet.settlement;

%% 2. Quarterly reset dates T_0, ..., T_16
N = 16;
resetDates    = zeros(1, N+1);
resetDates(1) = t0;
for i = 1:N
    resetDates(i+1) = addtodate(t0, 3*i, 'month');
end

delta    = yearfrac(resetDates(1:end-1), resetDates(2:end), 3); % Act/365, size 16
discGrid = linearRateInterp(dates, discounts, t0, resetDates);  % size 17
L0       = (discGrid(1:end-1) ./ discGrid(2:end) - 1) ./ delta; % size 16

%% 3. ATM strikes: forward swap rate for each cap maturity
maturitiesYears = [1, 2, 3, 4];
K_ATM = zeros(1, 4);
for m = 1:4
    nQ      = 4 * m;
    B_Tm    = discGrid(nQ + 1);
    BPV     = sum(delta(1:nQ) .* discGrid(2:nQ+1));
    K_ATM(m) = (1 - B_Tm) / BPV;
end

%% 4. Read flat vols from market table and interpolate at ATM strikes
strikes_table = [1.50,1.75,2.00,2.25,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,10.0] / 100;
vols_table = [
    0.140,0.130,0.129,0.121,0.133,0.138,0.144,0.150,0.172,0.191,0.202,0.216,0.239;
    0.224,0.197,0.175,0.180,0.192,0.204,0.210,0.214,0.223,0.236,0.249,0.261,0.281;
    0.238,0.217,0.200,0.198,0.203,0.205,0.208,0.214,0.229,0.243,0.256,0.267,0.282;
    0.242,0.224,0.209,0.204,0.204,0.202,0.202,0.205,0.217,0.229,0.240,0.250,0.266;
];

flatVolsATM_annual = zeros(1, 4);
for m = 1:4
    flatVolsATM_annual(m) = interp1(strikes_table, vols_table(m,:), K_ATM(m), ...
                                    'linear', 'extrap');
end

%% 5. Bootstrap spot vols nu_i from flat vols (linear constraint between annual maturities)
% Between two annual maturities, nu_i is linearly interpolated.
% One unknown per annual segment: nu at the right boundary (nu_left fixed from previous).

nu      = zeros(1, 15);
nu_left = interp1(maturitiesYears, flatVolsATM_annual, ...
                  yearfrac(t0, resetDates(2), 3), 'linear', 'extrap');

for k = 1:4
    i_left  = (k-1)*4 + 1;
    i_right = min(k*4, 15);

    mat_right = yearfrac(t0, resetDates(i_right+1), 3);
    mat_left  = yearfrac(t0, resetDates(i_left), 3);
    fv_right  = interp1(maturitiesYears, flatVolsATM_annual, mat_right, 'linear','extrap');
    fv_left   = interp1(maturitiesYears, flatVolsATM_annual, mat_left,  'linear','extrap');

    price_right = capPriceBMM_flat(fv_right, discGrid, delta, L0, resetDates, t0, i_right);
    if k == 1
        price_target = price_right;
    else
        price_left   = capPriceBMM_flat(fv_left, discGrid, delta, L0, resetDates, t0, i_left-1);
        price_target = price_right - price_left;
    end

    obj          = @(nu_r) segmentPrice(nu_r, nu_left, i_left, i_right, ...
                                        discGrid, delta, L0, resetDates, t0) - price_target;
    nu_right_sol = fzero(obj, fv_right);

    for i = i_left:i_right
        nu(i) = nu_left + (i - i_left) / (i_right - i_left + 1) * (nu_right_sol - nu_left);
    end
    nu_left = nu_right_sol;
end

%% 6. Correlation matrix rho_{ij} = exp(-lambda * |T_i - T_j|)
lambda = 0.1;
tGrid  = yearfrac(t0, resetDates(2:end), 3); % T_1..T_16, size 16
RHO    = exp(-lambda * abs(tGrid' - tGrid));  % 16x16
C      = chol(RHO, 'lower');

%% 7. Monte Carlo under spot measure
% B_i(t) is lognormal; chain advances reset date by reset date.
% Stochastic discount: D(t0, T_{k+1}) = prod_{j=0}^{k} B_j(T_j) accumulated along path.

Nsim    = 100000;
payoffs = zeros(Nsim, 1);

for sim = 1:Nsim

    % Initial forward bond values: B_i(t0) = disc(T_i)/disc(T_{i+1})
    B_sim = discGrid(1:end-1) ./ discGrid(2:end); % size 16

    % Stochastic discount factor accumulated along the path
    stochDisc = 1.0;

    for k = 0:N-1
        dt    = delta(k+1);
        Z_cor = C * randn(N, 1);

        % Update all bonds not yet fixed (i >= k+2)
        for i = k+2:N
            nu_i = nu(min(i-1, 15));

            % Convexity drift under spot measure: -sum_{j=k+1}^{i-1} rho_{ij}*nu_j*nu_i*dt
            drift_i = 0;
            for j = k+2:i-1
                nu_j    = nu(min(j-1, 15));
                drift_i = drift_i - RHO(i,j) * nu_j * nu_i * dt;
            end

            B_sim(i) = B_sim(i) * exp(drift_i - 0.5*nu_i^2*dt + nu_i*sqrt(dt)*Z_cor(i));
        end

        % Accumulate stochastic discount: multiply by B_{k+1}(T_k) just fixed
        stochDisc = stochDisc * B_sim(k+1);

        % Exotic caplet payoff: fixing at T_k, payment at T_{k+2}
        % First payment at 6M => k >= 1, last fixing at T_15 => k <= 14
        if k >= 1 && k <= 14
            i_curr = k + 1;
            i_prev = k;

            L_curr = (1/B_sim(i_curr) - 1) / delta(i_curr);
            L_prev = (1/B_sim(i_prev) - 1) / delta(i_prev);

            payoff_k = delta(i_curr) * max(L_curr - L_prev - 0.0005, 0);

            % Discount to t0: payment at T_{k+2} => divide by one more bond
            disc_to_payment = stochDisc * B_sim(i_curr+1);
            payoffs(sim)    = payoffs(sim) + payoff_k / disc_to_payment;
        end
    end
end

%% 8. Results
price    = mean(payoffs);
std_err  = std(payoffs) / sqrt(Nsim);

fprintf('Exotic Cap price (BMM): %.6f\n', price);
fprintf('Standard error:         %.6f\n', std_err);
fprintf('95%% CI: [%.6f, %.6f]\n', price - 1.96*std_err, price + 1.96*std_err);




