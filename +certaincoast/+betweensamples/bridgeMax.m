function M = bridgeMax(Y, sigma2, dt, U)
%BRIDGEMAX Maximum of each path including the exact Brownian-bridge maximum inside every time step.
%   M = certaincoast.betweensamples.bridgeMax(Y, SIGMA2, DT) returns, for each row of Y (paths sampled every DT), the
%   maximum over the whole time interval: inside each step the path is a Brownian bridge with local variance SIGMA2
%   per unit time, whose maximum is drawn exactly. Crossings that happen between grid points are therefore counted.
%   M = certaincoast.betweensamples.bridgeMax(Y, SIGMA2, DT, U) uses the given uniform numbers (size of Y without its
%   last column).
%
%   For a standardised Ornstein-Uhlenbeck path with rate KAPPA, SIGMA2 = 2 * KAPPA.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

a = Y(:, 1:end - 1);
b = Y(:, 2:end);
if nargin < 4 || isempty(U)
    U = rand(size(a));
end
m = 0.5 * (a + b + sqrt((b - a) .^ 2 - 2 * sigma2 * dt * log(U)));
M = max(m, [], 2);
end
