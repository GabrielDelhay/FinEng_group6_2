function fig = comparison_plot_ex4(x_vec, C_quad, C_fft1, C_fft2, C_fft3, C_mc)
% comparison_plot_ex4 - Plot comparison of Lewis pricing methods for Ex4

fig = figure('Position', [100 100 800 500]);
plot(x_vec*100, C_quad , 'k-' , 'LineWidth', 2 , 'DisplayName', 'Quadrature'); hold on;
plot(x_vec*100, C_fft1 , 'bo' , 'MarkerSize', 6 , 'DisplayName', 'FFT  dz = 0.0025');
plot(x_vec*100, C_fft2 , 'r+' , 'MarkerSize', 6 , 'DisplayName', 'FFT  freq extr = 500');
plot(x_vec*100, C_fft3 , 'g.' , 'MarkerSize', 12, 'DisplayName', 'FFT  dz = 0.01');
plot(x_vec*100, C_mc(:), 'm--', 'LineWidth', 1 , 'DisplayName', 'Monte Carlo');

xlabel('Moneyness x = log(F_0/K) [%]');
ylabel('Call price');
title('Ex4 - NIG: Lewis pricing methods');
legend('Location', 'northwest');
grid on;

end