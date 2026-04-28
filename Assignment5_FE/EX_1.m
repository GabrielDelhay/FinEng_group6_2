%% 1. Case Study: Structured bond

% Set strikes
strikes = [1.50, 1.75, 2.00, 2.25, 2.50, 3.00, 3.50, 4.00, 5.00, 6.00, 7.00, 8.00, 10.00];

% Set maturities (years)
maturities = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20];
n_strikes   = length(strikes);
n_maturities = length(maturities);


% Init volatility matrix
implied_vols = [
    14.0, 13.0, 12.9, 12.1, 13.3, 13.8, 14.4, 15.0, 17.2, 19.1, 20.2, 21.6, 23.9;
    22.4, 19.7, 17.5, 18.0, 19.2, 20.4, 21.0, 21.4, 22.3, 23.6, 24.9, 26.1, 28.1;
    23.8, 21.7, 20.0, 19.8, 20.3, 20.5, 20.8, 21.4, 22.9, 24.3, 25.6, 26.7, 28.2;
    24.2, 22.4, 20.9, 20.4, 20.4, 20.2, 20.2, 20.5, 21.7, 22.9, 24.0, 25.0, 26.6;
    24.3, 22.6, 21.2, 20.6, 20.4, 19.8, 19.5, 19.6, 20.5, 21.5, 22.6, 23.5, 25.0;
    24.3, 22.7, 21.4, 20.7, 20.2, 19.4, 18.9, 18.8, 19.3, 20.2, 21.2, 22.0, 23.5;
    24.1, 22.6, 21.4, 20.7, 20.1, 19.1, 18.4, 18.1, 18.4, 19.1, 20.0, 20.8, 22.2;
    23.9, 22.5, 21.4, 20.6, 20.0, 18.8, 18.0, 17.6, 17.6, 18.2, 19.0, 19.8, 21.1;
    23.7, 22.4, 21.3, 20.5, 19.8, 18.5, 17.6, 17.1, 17.0, 17.6, 18.3, 19.0, 20.3;
    23.5, 22.2, 21.2, 20.4, 19.6, 18.3, 17.3, 16.8, 16.5, 16.9, 17.6, 18.3, 19.5;
    23.0, 21.7, 20.8, 20.0, 19.3, 17.9, 16.9, 16.2, 15.8, 16.0, 16.5, 17.1, 18.1;
    22.3, 21.2, 20.3, 19.5, 18.7, 17.3, 16.3, 15.5, 15.0, 15.1, 15.5, 16.0, 16.9;
    21.6, 20.4, 19.5, 18.8, 18.0, 16.6, 15.5, 14.7, 14.1, 14.1, 14.5, 15.0, 15.9
];

implied_vols = implied_vols / 100;   % -> convert vols to decimals
% Load and bootstrap the zero rate curve
formatDate = 'dd/mm/yyyy';
[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap.xls', formatDate);
[dates, discounts, zeroRates] = bootstrap(datesSet, ratesSet);
t0 = dates(1);  % Settlement date

% Exercise 1a: LMM Spot Vol Calibration
%  
%  GOAL: Extract LMM spot vols from market flat cap vols
%
%  THEORY:
%    - A Cap(T,K) = sum of Caplets, each priced with Black:
%        caplet_i = B(0,T_{i+1}) * delta_i * [L_i*N(d1) - K*N(d2)]
%      where d1,d2 use the SPOT vol sigma_i (not flat vol)
%    - Flat vol Sigma_n: same vol for ALL caplets in cap -> 
%        Cap(T_n, K) = sum_i caplet(i, Sigma_n)   [market quote]
%    - Spot vol sigma_i: TRUE vol of forward rate i ->
%        Cap(T_n, K) = sum_i caplet(i, sigma_i)   [model decomp]
%    - Bootstrap: for each new maturity T_n, the price of the
%      NEW caplets = Cap(T_n) - Cap(T_{n-1}), so we solve for
%      the spot vols in the new bucket via the linear constraint
%      (slide 23):
%        sigma_i = sigma(T_alpha) + (T_i - T_alpha)/(T_beta - T_alpha)
%                  * (sigma(T_beta) - sigma(T_alpha))
%      i.e. spot vols are LINEARLY interpolated within each bucket
%
% the caplet from T_0 to T_1 (first 3m) is assumed already fixed at t_0, so caps start from caplet 2.

%%  SECTION 1: Build the quarterly caplet schedule
%  The cap schedule is quarterly from t0.
%  dates(1) = t0 = settlement
%  dates(2) = t0 + ~3m (first quarterly date)
%  We select from the bootstrapped dates only the quarterly
%  ones that fall within the cap maturities (up to 20y).

% The bootstrapped dates vector already contains quarterly dates.
% We use ALL dates as potential caplet dates.
% T_i are the reset dates, T_{i+1} are the payment dates.

% Last cap maturity is 20y from t0
max_maturity_date = t0 + 20*365 + 5; % approximate upper bound

% Select dates within range (exclude t0 itself)
cap_dates = dates(dates > t0 & dates <= max_maturity_date);
n_cap_dates = length(cap_dates);

% The bootstrapped dates are NOT purely quarterly. We therefore
% interpolate discount factors from the bootstrap curve onto
% the quarterly caplet schedule defined above.
% log-linear interpolation is used (standard for discount factors).

%% SECTION 2: Compute discount factors and forward rates on the caplet schedule
% Interpolate discount factors at caplet dates
% We interpolate log-linearly (standard for discount factors)
log_discounts = log(discounts);
B_cap = exp(interp1(dates, log_discounts, cap_dates, 'linear'));
all_cap_dates    = [t0; cap_dates(:)];
all_B            = [1;  B_cap(:)];   %B(0, t0) = 1 by definition

n_caplets = length(cap_dates); % number of potential payment dates

% Pre-allocate forward rates and day count fractions
fwd_rates = zeros(n_caplets, 1); % L_i = forward rate for caplet i
delta_fwd = zeros(n_caplets, 1); % Act/360 year fraction for coupon
tau_expiry = zeros(n_caplets, 1); % Act/365 year fraction for Black

% I DAYCOUNT DEVO RIGUARDARLI

for i = 1:n_caplets
    T_reset   = all_cap_dates(i);    % = T_i
    T_payment = all_cap_dates(i+1);  % = T_{i+1}
    
    % !! DAYCOUNT: Act/360 for delta_i (coupon accrual)
    delta_fwd(i) = (T_payment - T_reset) / 360;
    
    % !! DAYCOUNT: Act/365 for T_i in Black vol scaling
    tau_expiry(i) = (T_reset - t0) / 365;
    
    % Forward Euribor 3m rate from discount factors
    % L_i = (1/delta_i) * (B(0,T_i)/B(0,T_{i+1}) - 1)
    B_Ti       = all_B(i);    % B(0, T_i)
    B_Ti1      = all_B(i+1);  % B(0, T_{i+1})
    fwd_rates(i) = (1/delta_fwd(i)) * (B_Ti/B_Ti1 - 1);
end

%% SECTION 3: Map cap maturities to caplet indices
%
%  Each cap maturity (1y, 2y, ..., 20y) corresponds to a set of caplets 
%  ending at that maturity. We find which caplet index corresponds to each 
%  maturity. (FIRST CAPLET EXCLUDED)
%  So Cap(T_n) starts from caplet index 2 (i=2).

% Find the caplet index for each cap maturity
% cap_maturity_idx(j) = last caplet index for cap j
cap_maturity_idx = zeros(n_maturities, 1);

for j = 1:n_maturities
    % Target payment date: t0 + maturity_j years (approx)
    target_date = t0 + maturities(j) * 365.25;
    
    % Find the caplet whose payment date is closest to target
    [~, idx] = min(abs(cap_dates - target_date));
    cap_maturity_idx(j) = idx;
end

i_start = 2;  % First caplet index

%% SECTION 4: Bootstrap LMM spot vols
%  ALGORITHM:
%  For each strike k:
%    For each maturity bucket [T_{alpha}, T_{beta}]:
%      1. Compute market cap price Cap_mkt(T_beta, K) using flat vol
%      2. Subtract already-priced caplets from previous buckets
%         -> get Delta_C = price of NEW caplets in this bucket
%      3. Within the bucket, spot vols are LINEARLY interpolated:
%           sigma_i = sigma_alpha + (T_i - T_alpha)/(T_beta - T_alpha)
%                     * (sigma_beta - sigma_alpha)
%         where sigma_alpha is known (last bucket's endpoint)
%         and sigma_beta is the unknown to solve for
%      4. Solve for sigma_beta such that 
%           sum_{i in bucket} caplet(i, sigma_i) = Delta_C
%         using fzero (1D root finding)
%
%  OUTPUT: spot_vols(i, k) = spot vol of caplet i for strike k
%          (same grid as flat vols: 13 strikes x n_caplets)

% Pre-allocate spot vol matrix
% Rows = caplet index (from i_start to last cap maturity)
% Cols = strike index
n_total_caplets = cap_maturity_idx(end); % last caplet index
spot_vols = zeros(n_total_caplets, n_strikes);

% Loop over strikes
for k = 1:n_strikes
    K = strikes(k) / 100; % convert to decimal
    
    % sigma at the left edge of the first bucket:
    % Before the first cap maturity, we have no information.
    % We initialize sigma_alpha = 0 (no vol before first cap).
    % The first bucket goes from i_start to cap_maturity_idx(1).
    sigma_alpha = 0; % left edge of first bucket (unknown convention,
                     % set to 0: effectively means first bucket solved DA RIGUARDARE
                     % with flat spot vol = sigma_beta throughout)
    
    % Running total of already-priced caplet prices (for subtraction)
    cap_price_already_priced = 0;
    
    % Index of the last processed caplet
    i_prev_end = i_start - 1;
    
    % Loop over maturity buckets
    for j = 1:n_maturities
        
        % Step 1: Market cap price at T_beta
        % Using flat vol Sigma_j for ALL caplets from i_start to end
        % Cap price = sum_{i=i_start}^{i_beta} caplet_Black(i, Sigma_j)
        Sigma_j = implied_vols(j, k); % flat vol for this maturity/strike
        i_beta  = cap_maturity_idx(j);    % last caplet index for this maturity
        
        cap_price_flat = 0;
        for i = i_start:i_beta
            cap_price_flat = cap_price_flat + ...
                caplet_black_LMM(fwd_rates(i), K, delta_fwd(i), ...
                                  all_B(i+1), tau_expiry(i), Sigma_j);
        end
        
        % Step 2: Delta cap = new caplets only
        % These are the caplets in (i_prev_end, i_beta]
        delta_cap_price = cap_price_flat - cap_price_already_priced;
        
        % Caplet indices in this new bucket
        i_alpha_new = i_prev_end + 1; % first new caplet
        i_beta_new  = i_beta;          % last new caplet
        
        % Step 3: Solve for sigma_beta 
        % Within the bucket, spot vols are linear in caplet index:
        %   sigma_i = sigma_alpha + (i-i_alpha)/(i_beta-i_alpha) 
        %             * (sigma_beta - sigma_alpha)
        T_alpha_new = tau_expiry(i_alpha_new); % time of first new caplet
        T_beta_new  = tau_expiry(i_beta_new);  % time of last new caplet
        
        % Define the function whose root we seek:
        % f(sigma_beta) = sum of caplets with linearly interpolated vols
        %                 - delta_cap_price = 0
        f = @(sigma_beta) ...
            sum_caplets_linear_vol(i_alpha_new, i_beta_new, ...
                fwd_rates, K, delta_fwd, all_B, tau_expiry, ...
                sigma_alpha, sigma_beta, T_alpha_new, T_beta_new) ...
            - delta_cap_price;
        
        % Solve with fzero (bracketing: vol must be positive)
        % Initial guess: sigma_alpha (or small positive value)
        sigma_beta_init = max(sigma_alpha, 0.001);
        
        try
            sigma_beta = fzero(f, [1e-6, 2.0]); % vol between 0 and 200%
        catch
            % If fzero fails (e.g. no sign change), use last known vol
            % This can happen for very short or very long maturities
            warning('fzero failed for maturity %dy, strike %.2f%%. Using sigma_alpha.', ...
                    maturities(j), strikes(k));
            sigma_beta = sigma_alpha;
        end
        
        % Step 4: Store spot vols for this bucket 
        for i = i_alpha_new:i_beta_new
            % Linear interpolation of spot vol within bucket 
            if T_beta_new > T_alpha_new
                spot_vols(i, k) = sigma_alpha + ...
                    (tau_expiry(i) - T_alpha_new) / (T_beta_new - T_alpha_new) ...
                    * (sigma_beta - sigma_alpha);
            else
                % Single caplet bucket: flat vol
                spot_vols(i, k) = sigma_beta;
            end
        end
        
        % Update for next iteration
        % The new sigma_alpha is sigma_beta of this bucket
        sigma_alpha = sigma_beta;
        
        % Update already-priced cap price:
        % Recompute with the NOW-KNOWN spot vols
        cap_price_already_priced = 0;
        for i = i_start:i_beta_new
            cap_price_already_priced = cap_price_already_priced + ...
                caplet_black_LMM(fwd_rates(i), K, delta_fwd(i), ...
                                  all_B(i+1), tau_expiry(i), spot_vols(i,k));
        end
        
        i_prev_end = i_beta_new;
    end % end maturity loop
    
end % end strike loop

%% SECTION 5: Display results
%  Show spot vols at cap maturity dates (same grid as input) for easy comparison with flat vols

fprintf('\n=== LMM Spot Vols (at cap maturity dates) [%%] ===\n');
fprintf('Maturity |');
for k = 1:n_strikes
    fprintf(' K=%.2f%%', strikes(k));
end
fprintf('\n');

for j = 1:n_maturities
    fprintf('  %3dy   |', maturities(j));
    i_beta = cap_maturity_idx(j);
    for k = 1:n_strikes
        fprintf('  %6.2f ', spot_vols(i_beta, k) * 100);
    end
    fprintf('\n');
end