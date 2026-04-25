function fig = plot_smile(strikes, smiles_mkt, F0)
% plot_smile - Plot the implied volatility smile with spline interpolation

% Extended fine grid for smooth smile plot
K_fine = linspace(strikes(1), strikes(end), 500);
vol_fine = interp1(strikes, smiles_mkt, K_fine, 'spline');

fig = figure('Name', 'Implied Volatility Smile - Eurostoxx 50 (15 Feb 2008)');
plot(strikes, smiles_mkt*100, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6); hold on;
plot(K_fine, vol_fine*100, 'b-', 'LineWidth', 1.5);
xline(F0, 'r--', 'LineWidth', 1.5, 'Label', sprintf('ATM F_0 = %.1f', F0));

ylabel('Implied Volatility (%)');
xlabel('Strike K');
title('Eurostoxx 50 Implied Volatility Smile - 1Y, 15 Feb 2008');
legend('Market vols', 'Spline interpolation', 'ATM price F0', 'Location', 'NorthEast');
grid on;

end