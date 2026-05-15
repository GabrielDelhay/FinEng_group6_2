function printPricingResults(priceSwaption, B0_10, B0_10_tree, upperBound, lowerBound)
% PRINTPRICINGRESULTS Prints the Bermudan swaption pricing results and checks theoretical bounds.
%
% Inputs:
%   priceSwaption - Price of the Bermudan payer swaption
%   B0_10         - Bond price B(0,T_N) from the initial yield curve
%   B0_10_tree    - Bond price B(0,T_N) computed from the pricing tree
%   upperBound    - Upper bound price (e.g., corresponding Cap)
%   lowerBound    - Lower bound price (e.g., max European swaption via Jamshidian)

fprintf('\n=== Exercise 2 ===\n');

% Print detailed results
fprintf('Bermudan payer swaption price : %.6f\n', priceSwaption);
fprintf('B(0,T_N) tree vs curve        : %.2e\n', B0_10 - B0_10_tree);
fprintf('Upper bound (cap)             : %.6f\n', upperBound);
fprintf('Lower bound (max Eur Jamsh.)  : %.6f\n', lowerBound);

% Print summary and check if the price is within bounds
fprintf('\n=== Summary ===\n');
fprintf('Lower    : %.6f\n', lowerBound);
fprintf('Bermudan : %.6f\n', priceSwaption);
fprintf('Upper    : %.6f\n', upperBound);

% Evaluate boolean condition and print as string
isInBand = priceSwaption >= lowerBound && priceSwaption <= upperBound;
fprintf('In band? : %s\n', string(isInBand));

end