function [sigma, a, K, N, nStepsPerYear, dt, dx, muHat, lMax, x] = getTreeParameters()
% GETTREEPARAMETERS Initializes and returns the model parameters and tree settings.

% Model parameters
sigma = 0.8e-2;
a     = 0.11;
K     = 0.05;
N     = 10;

% Tree settings
nStepsPerYear = 365;
dt    = 1/nStepsPerYear;
dx    = sigma * sqrt(3*dt);
muHat = 1 - exp(-a*dt);
lMax  = ceil((1 - sqrt(2/3))/muHat);
x     = dx * (lMax:-1:-lMax)';

end