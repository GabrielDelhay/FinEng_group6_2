function [err1, err2, err3, errMC] = errors_ex4(C_quad, C_fft1, C_fft2, C_fft3, C_mc)
% errors_ex4 - Compute and print max absolute errors vs Quadrature benchmark

err1 = max(abs(C_quad - C_fft1));
err2 = max(abs(C_quad - C_fft2));
err3 = max(abs(C_quad - C_fft3));
errMC = max(abs(C_quad - C_mc(:)));

fprintf('Max FFT error vs Quadrature:\n');
fprintf('   dz = 0.0025      : %.3e\n', err1);
fprintf('   freq extr = 500  : %.3e\n', err2);
fprintf('   dz = 0.01        : %.3e\n', err3);
fprintf('Max MC error vs Quadrature: %.3e\n', errMC);

end