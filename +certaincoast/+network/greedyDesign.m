function [S, gain, steps] = greedyDesign(Gcand, H, GamF, noiseVar, budget, cost)
%GREEDYDESIGN Cost-aware greedy choice of sensors for the largest goal-oriented information gain.
%   [S, GAIN] = certaincoast.network.greedyDesign(GCAND, H, GAMF, NOISEVAR, BUDGET) adds, one at a time, the
%   candidate (row of GCAND) with the largest gain in eigOperator per unit cost, until the budget is spent. With no
%   COST every candidate costs 1, so BUDGET is the number of sensors. NOISEVAR is a scalar or one variance per
%   candidate. S lists the chosen rows in the order they were added, GAIN is the information gain of the design.
%   [S, GAIN, STEPS] = certaincoast.network.greedyDesign(..., COST) uses the given costs and also returns, for every
%   step, the candidate added, the gain of the design and the marginal gain.
%
%   The goal-oriented objective is only weakly submodular, so the classical (1 - 1/e) guarantee does not hold in
%   general; compare with an exhaustive search for small budgets (poolingLoss does).

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

M = size(Gcand, 1);
if nargin < 6 || isempty(cost)
    cost = ones(M, 1);
end
if isscalar(noiseVar)
    noiseVar = noiseVar * ones(M, 1);
end
cost = cost(:);
noiseVar = noiseVar(:);
S = [];
spent = 0;
gain = 0;
steps = struct('added', {}, 'gain', {}, 'marginal', {});
remaining = true(M, 1);
while true
    bestJ = 0;
    bestRatio = -Inf;
    bestGain = gain;
    for j = find(remaining).'
        if spent + cost(j) > budget
            continue
        end
        T = [S, j];
        e = certaincoast.network.eigOperator(Gcand(T, :), H, GamF, diag(noiseVar(T)));
        if (e - gain) / cost(j) > bestRatio
            bestRatio = (e - gain) / cost(j);
            bestJ = j;
            bestGain = e;
        end
    end
    if bestJ == 0
        break
    end
    S = [S, bestJ]; %#ok<AGROW>
    spent = spent + cost(bestJ);
    steps(end + 1) = struct('added', bestJ, 'gain', bestGain, 'marginal', bestGain - gain); %#ok<AGROW>
    remaining(bestJ) = false;
    gain = bestGain;
end
end
