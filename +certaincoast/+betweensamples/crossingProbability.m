function P = crossingProbability(kappa, z, horizon, dt, nPaths)
%CROSSINGPROBABILITY Probability that the indicator crosses a limit at any time in a season, not only at the samples.
%   P = certaincoast.betweensamples.crossingProbability(KAPPA, Z, HORIZON) returns a numel(KAPPA)-by-numel(Z) array:
%   the probability that a stationary standardised Ornstein-Uhlenbeck indicator with rate KAPPA (1/days) reaches the
%   standardised limit Z at some time in [0, HORIZON] (days). For chloride with mean MU and standard deviation SD and
%   a limit L, Z = (L - MU) / SD.
%   P = certaincoast.betweensamples.crossingProbability(KAPPA, Z, HORIZON, DT, NPATHS) sets the simulation step
%   (default 0.25 day) and the number of paths (default 4000).
%
%   The same random numbers are used for every KAPPA (common random numbers), so the curve over KAPPA is smooth. The
%   probability does not decrease with KAPPA (Slepian's inequality): the marginal law of the samples alone gives only
%   a lower bound, 1 - Phi(Z), and no upper bound (Proposition S1 of the preparatory experiment PATH-GO).
%
%   See also certaincoast.betweensamples.timeScaleInterval, certaincoast.betweensamples.crossingBand.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 4 || isempty(dt)
    dt = 0.25;
end
if nargin < 5 || isempty(nPaths)
    nPaths = 4000;
end
n = round(horizon / dt);
E = randn(nPaths, n + 1);
U = rand(nPaths, n);
z = z(:).';
P = zeros(numel(kappa), numel(z));
for i = 1:numel(kappa)
    X = certaincoast.betweensamples.ouPaths(kappa(i), nPaths, n, dt, E);
    M = certaincoast.betweensamples.bridgeMax(X, 2 * kappa(i), dt, U);
    P(i, :) = mean(bsxfun(@ge, M, z), 1);
end
end
