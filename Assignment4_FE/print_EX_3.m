function print_EX_3( x_vec3, C_quad3, C_res3, C_mc3, C_fft3)
    fprintf('Moneyness:    %10.5f  %10.5f  %10.5f\n', x_vec3);
    fprintf('Quadrature:   %10.4f  %10.4f  %10.4f\n', C_quad3);
    fprintf('Residuals:    %10.4f  %10.4f  %10.4f\n', C_res3);
    fprintf('MonteCarlo:   %10.4f  %10.4f  %10.4f\n', C_mc3);
    fprintf('FFT:          %10.4f  %10.4f  %10.4f\n\n', C_fft3);
end
