function X = ouPaths(kappa, nPaths, nSteps, dt, epsilon)
%OUPATHS Stationary standardised Ornstein-Uhlenbeck paths with exact transitions.
%   X = certaincoast.betweensamples.ouPaths(KAPPA, NPATHS, NSTEPS, DT) returns an NPATHS-by-(NSTEPS + 1) array of paths
%   of dX = -KAPPA X dt + sqrt(2 KAPPA) dW (mean 0, variance 1, correlation time 1/KAPPA), started in the stationary
%   law and sampled every DT.
%   X = certaincoast.betweensamples.ouPaths(..., EPSILON) uses the given standard normal numbers
%   (NPATHS-by-(NSTEPS + 1)), so that several values of KAPPA can share the same random numbers.
%
%   See also certaincoast.betweensamples.bridgeMax, certaincoast.betweensamples.crossingProbability.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 5 || isempty(epsilon)
    epsilon = randn(nPaths, nSteps + 1);
end
rho = exp(-kappa * dt);
s = sqrt(1 - rho ^ 2);
X = zeros(size(epsilon));
X(:, 1) = epsilon(:, 1);
for k = 1:nSteps
    X(:, k + 1) = rho * X(:, k) + s * epsilon(:, k + 1);
end
end
