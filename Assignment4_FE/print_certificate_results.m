function print_certificate_results(floaterPV, spreadBPV, protectionPV, M, C_basket, SE_basket, CI_basket, X, SE_X, CI_X, N)
% print_certificate_results - Print deterministic components, basket option price, and upfront X

fprintf('Deterministic components (per unit of N):\n  floaterPV     = %.6f\n  spreadBPV     = %.6f\n  protectionPV  = %.6f\n\n', floaterPV, spreadBPV, protectionPV);

fprintf('Basket call price (MC, M = %.0e):\n  C_basket = %.6f   (SE = %.6f)\n  95%% CI   = [%.6f , %.6f]\n\n', M, C_basket, SE_basket, CI_basket(1), CI_basket(2));

fprintf('Fair upfront X%%:\n  X       = %.8f %%   (SE = %.6f %%)\n  95%% CI   = [%.6f %% , %.6f %%]\n', 100*X, 100*SE_X, 100*CI_X(1), 100*CI_X(2));

fprintf('  cash    = %.8f M EUR   (on N = %.0f M EUR)\n', N*X/1e6, N/1e6);

end