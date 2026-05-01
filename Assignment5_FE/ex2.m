%% Exotic Cap - Calibration & Pricing under Bond Market Model (BMM)

clc; close all; clear all;
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
    raw_date         = addtodate(t0, 3*i, 'month');
    resetDates(i+1)  = busdate(raw_date, 'modifiedfollow');
end

delta    = yearfrac(resetDates(1:end-1), resetDates(2:end), 3); % Act/365
B0_T = linearRateInterp(dates, discounts, t0, resetDates); 
L0       = (B0_T(1:end-1) ./ B0_T(2:end) - 1) ./ delta;

%% 3. ATM strikes: forward swap rate for each cap maturity
maturitiesYears = [1, 2, 3, 4];
K_ATM = zeros(1, 4);
for m = 1:4
    nQ      = 4 * m;
    B_Tm    = B0_T(nQ + 1);
    BPV     = sum(delta(1:nQ) .* B0_T(2:nQ+1));
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

%% Spot vol bootstrap (BMM, single strike K_i = L0(i) per caplet)
nu = spotvolbootstrap(maturitiesYears, flatVolsATM_annual, t0, resetDates, B0_T, delta, L0);

%% 6. Correlation matrix rho_{ij} = exp(-lambda * |T_i - T_j|)
lambda = 0.1;
tGrid  = yearfrac(t0, resetDates(2:end), 3); % T_1..T_16
RHO    = exp(-lambda * abs(tGrid' - tGrid));  % 16x16
C      = chol(RHO, 'lower');

%% 7. Monte Carlo under Spot Measure (Bond Market Model)

Nsim    = 100000;
payoffs = zeros(Nsim,1);


for sim = 1:Nsim
    % B_sim(i) = P(t0,T_i)/P(t0,T_{i-1}) < 1, i=1..16 (forward discount factor)
    B_sim = B0_T(2:end) ./ B0_T(1:end-1);
    fixedLibors = zeros(1, N);
    fixedLibors(1) = L0(1);  % L_1 known at t0 (deterministic)
    cumDiscount = 1.0;

    for k = 0:N-1
        dt    = delta(k+1);            % time step [T_k, T_{k+1}]
        Z     = randn(N, 1);
        dW    = sqrt(dt) * (C * Z);    % correlated Brownian increments (16 x 1)
        B_old = B_sim;

        % Evolve forward bonds B_{k+2},...,B_N from T_k to T_{k+1}.
        % Bond B_i = P(t,T_i)/P(t,T_{i-1}), vol nu(i-1), noise dW(i-1).
        for i = k+2:N
            nu_i  = nu(i-1);
            % Spot-measure convexity drift: -nu_i * sum_{j=k+2}^{i-1} rho_{ij}*nu_j
            drift = 0.0;
            for j = k+2:i-1
                nu_j  = nu(j-1);
                drift = drift - RHO(i-1, j-1) * nu_i * nu_j;
            end
            B_sim(i) = B_old(i) * exp((drift - 0.5*nu_i^2)*dt - nu_i*dW(i-1));
        end

        % Fix L_{k+2} at T_{k+1}: spot Libor for period [T_{k+1}, T_{k+2}]
        if k <= N-2
            fixedLibors(k+2) = (1/B_sim(k+2) - 1) / delta(k+2);
        end

        % Accumulate path discount factor: D(t0, T_{k+1}) = prod_{j=0}^{k} P(T_j, T_{j+1})
        cumDiscount = cumDiscount * B_old(k+1);

        % Exotic caplet payoff at T_{k+2}:
        %   delta_{k+1} * max( L_{k+1}(T_k) - L_k(T_{k-1}) - spread, 0 )
        % First payment at T_3 (k=1), last at T_16 (k=14).
        if k >= 1 && k <= N-2
            L_prev = fixedLibors(k);    % L_k fixed at T_{k-1}
            L_curr = fixedLibors(k+1);  % L_{k+1} fixed at T_k
            payoff = delta(k+1) * max(L_curr - L_prev - 0.0005, 0.0);
            % Discount from T_{k+2} to t0: cumDiscount*B_sim(k+2) = P(t0,T_{k+2})
            pv = cumDiscount * B_sim(k+2) * payoff;
            payoffs(sim) = payoffs(sim) + pv;
        end
    end
end

%% 8. Results
price   = mean(payoffs);
std_err = std(payoffs)/sqrt(Nsim);

fprintf('Exotic Cap price (BMM): %.8f\n', price);
fprintf('Standard error:         %.8f\n', std_err);
fprintf('95%% CI: [%.8f, %.8f]\n', price - 1.96*std_err, price + 1.96*std_err);