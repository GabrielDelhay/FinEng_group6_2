function fig = plot_calibration(F0, strikes, smiles_mkt, IV_model, RMSE)
% plot_calibration - Plot market vs model implied volatility for global calibration

fig = figure;
x_mkt = log(F0 ./ strikes);
plot(x_mkt, smiles_mkt*100, 'b o-', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', 'Market');
hold on;
plot(x_mkt, IV_model*100, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('nMV \\alpha=2/3 (RMSE=%.3f%%)', RMSE));

xlabel('Moneyness x = ln(F_0/K)');
ylabel('Implied Volatility (%)');
title('Global Calibration – nMV \alpha=2/3');
legend('Location', 'best');
grid on;

end