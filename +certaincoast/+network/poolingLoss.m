function out = poolingLoss(Gcand, H, blocks, cross, noiseVar, budget, rhoGrid)
%POOLINGLOSS What is lost by designing the network as if the facies were one: the poolability phase boundary.
%   OUT = certaincoast.network.poolingLoss(GCAND, H, BLOCKS, CROSS, NOISEVAR, BUDGET, RHOGRID) computes, for each
%   cross-facies correlation RHOF in RHOGRID, the loss
%
%       L(RHOF) = max_d EIG_RHOF(d) - EIG_RHOF(d0),
%
%   where the best design d of BUDGET sensors is found by exhaustive search under the true prior faciesPrior(RHOF)
%   and d0 is the best design under the pooled prior faciesPrior(1). OUT has fields rho, L (nats), Lambda = L / L(0),
%   rhoC (the largest RHOF with a loss above round-off: above it, pooling costs nothing) and d0.
%   Use [] for CROSS to take the default coupling of faciesPrior.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

combs = nchoosek(1:size(Gcand, 1), budget);
d0 = bestDesign(Gcand, H, certaincoast.network.faciesPrior(1, blocks, cross), noiseVar, combs);
L = zeros(numel(rhoGrid), 1);
for i = 1:numel(rhoGrid)
    Gam = certaincoast.network.faciesPrior(rhoGrid(i), blocks, cross);
    dStar = bestDesign(Gcand, H, Gam, noiseVar, combs);
    L(i) = max(designGain(Gcand, H, Gam, noiseVar, dStar) - designGain(Gcand, H, Gam, noiseVar, d0), 0);
end
Lambda = L / max(L(1), eps);
active = L > 1e-12 * max(L(1), eps);
rhoC = NaN;
if any(active)
    rhoC = rhoGrid(find(active, 1, 'last'));
end
out = struct('rho', rhoGrid(:), 'L', L, 'Lambda', Lambda, 'rhoC', rhoC, 'd0', d0);
end

function d = bestDesign(Gcand, H, Gam, noiseVar, combs)
best = -Inf;
d = combs(1, :);
for c = 1:size(combs, 1)
    e = designGain(Gcand, H, Gam, noiseVar, combs(c, :));
    if e > best
        best = e;
        d = combs(c, :);
    end
end
end

function e = designGain(Gcand, H, Gam, noiseVar, S)
if isscalar(noiseVar)
    Gn = noiseVar * eye(numel(S));
else
    Gn = diag(noiseVar(S));
end
e = certaincoast.network.eigOperator(Gcand(S, :), H, Gam, Gn);
end
