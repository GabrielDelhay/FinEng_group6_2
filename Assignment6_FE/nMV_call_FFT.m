function C = nMV_call_FFT(p, alpha, F0, K_vec, B, T)
    % Lewis via FFT
    %
    % Integral:  I(x) = int_{-inf}^{+inf} dxi/(2pi) * e^{-i*xi*x}
    %                   * phi(-xi-i/2) / (xi^2 + 1/4)
    %
    % Discretization:
    %   xi_j = xi_1 + (j-1)*dxi,   j = 1,...,N
    %   x_k  = x_1  + (k-1)*dx
    %   dxi * dx = 2*pi/N
    
    M   = 14;
    N   = 2^M;
    dxi = 0.01;         % frequency step
    dx  = 2*pi/(N*dxi); % moneyness step
    
    % xi grid: centered around 0  (xi goes from negative to positive)
    j    = (0:N-1)';
    xi1  = -(N/2)*dxi;           % starting frequency
    xi_j = xi1 + j*dxi;         % xi_j = xi_1 + (j-1)*dxi
    
    % x grid: centered around 0
    x1   = -(N/2)*dx;
    x_k  = x1 + j*dx;
    
    % Integrand evaluated on xi grid (this is what we FFT-transform)
    % phi must be analytic in strip Im(xi) in [-1,0]: OK by construction
    phi_vals  = char_fun(-xi_j - 1i/2, p, alpha, T);   % phi(-xi-i/2)
    integrand = phi_vals ./ (xi_j.^2 + 0.25);           % divide by xi^2+1/4
    
    % To use FFT we write
    %   I(x_k) = (dxi/2pi) * sum_j integrand(xi_j) * e^{-i*xi_j*x_k}
    %          = (dxi/2pi) * e^{-i*xi_1*x_k} * FFT[ e^{-i*(j-1)*dxi*x_1} * integrand_j ]
    %
    % Define f_j (the sequence fed into FFT):
    f_j = integrand .* exp(-1i * xi_j * x1);  % multiply by e^{-i*xi_j*x_1}
    % Note: xi_j*x1 = (xi_1 + (j-1)*dxi)*x1
    %               = xi_1*x1  +  (j-1)*dxi*x1
    % so f_j = integrand * e^{-i*xi_1*x1} * e^{-i*(j-1)*dxi*x1}
    % the e^{-i*xi_1*x1} is a global phase (const), absorbed in prefactor
    
    % FFT
    FFT_out = fft(f_j);   % FFT_out(k) = sum_j f_j * e^{-2pi*i*(j-1)*(k-1)/N}
    
    % Prefactor to convert FFT output to the integral at x_k 
    % I(x_k) = (dxi/2pi) * e^{-i*xi_1*x_k} * FFT_out(k)
    prefactor = exp(-1i * xi1 * x_k);
    I_k = (dxi / (2*pi)) * prefactor .* FFT_out;
    
    % Lewis formula: c(x)/(B*F0) = 1 - e^{-x/2} * Re{ I(x) }
    C_grid = B * F0 * (1 - exp(-x_k/2) .* real(I_k));
    C_grid = max(C_grid, 0);
    
    % Interpolate at requested strikes
    x_req = log(F0 ./ K_vec(:));
    
    % Check all requested moneyness values fall inside the grid
    if any(x_req < x_k(1)) || any(x_req > x_k(end))
        warning('Some strikes fall outside FFT moneyness grid! Increase M or adjust dxi.');
    end
    
    C = interp1(x_k, C_grid, x_req, 'spline');
    C = max(real(C), 0);
end
