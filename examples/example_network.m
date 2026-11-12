%EXAMPLE_NETWORK Where to put three salinity sensors in a coastal section: for the salt mass, or for a supply well?
%   Synthetic Henry section (16 x 8 cells) with an uncertain log-conductivity field (standard deviation 0.5,
%   correlation length 0.4 m). The sensitivities of 12 candidate salinity sensors, of the salt mass near the toe and
%   of the salinity at a supply cell come from two adjoint solves. The example compares the greedy design with the
%   exhaustive optimum, and the best design for the salt mass with the best design for the question "is the supply
%   cell above its limit?". Synthetic data only. Run from the repository root (a few seconds).

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

addpath(fileparts(fileparts(mfilename('fullpath'))));
g = certaincoast.model.henryGrid(16, 8);
P = certaincoast.model.henryParameters();
theta = zeros(g.N, 1);
U = certaincoast.model.solveDensity(theta, g, P, 1e-12);
[J, ~, Jtheta] = certaincoast.model.densityJacobian(U, theta, g, P);

[X, Y] = ndgrid(g.xc, g.yc);                    % cell centres in the flat order of the state
X = X(:);
Y = Y(:);
cellAt = @(x, y) find(abs(X - x) < g.dx / 2 & abs(Y - y) < g.dy / 2);
cand = zeros(12, 1);
k = 0;
for x = [0.5625 1.0625 1.4375 1.8125]
    for y = [0.0625 0.4375 0.8125]
        k = k + 1;
        cand(k) = cellAt(x, y);
    end
end
supply = cellAt(1.5625, 0.3125);                % a mid-depth well 0.44 m from the sea
zone = find(X > 1.0 & X < 1.6 & Y < 0.25);       % the toe zone

dOdu = sparse(1:12, g.N + cand, 1, 12, 2 * g.N);
dqdu = [zeros(1, g.N), zeros(1, g.N)];
dqdu(g.N + zone) = g.dx * g.dy;                  % salt mass in the toe zone
dsdu = sparse(1, g.N + supply, 1, 1, 2 * g.N);   % salinity at the supply cell
[G, H] = certaincoast.network.adjointSensitivities(J, Jtheta, [dOdu; dsdu], dqdu);
Hs = G(13, :);
G = G(1:12, :);

D = sqrt(bsxfun(@minus, X, X.') .^ 2 + bsxfun(@minus, Y, Y.') .^ 2);
Gam = 0.25 * exp(-D / 0.4);                      % prior covariance of the log-conductivity
noise = 0.02 ^ 2;                                % salinity sensor noise (2% of sea-water salinity)

[S, gain] = certaincoast.network.greedyDesign(G, H, Gam, noise, 3);
combs = nchoosek(1:12, 3);
gains = zeros(size(combs, 1), 1);
Sigma = [G; -Hs] * Gam * [G; -Hs].';             % joint covariance of the 12 sensors and of -salinity at supply
sdS = sqrt(Sigma(13, 13));
limit = U(g.N + supply) + 0.3 * sdS;             % the supply cell is "salinised" above this value
risks = zeros(size(combs, 1), 1);
for c = 1:size(combs, 1)
    gains(c) = certaincoast.network.eigOperator(G(combs(c, :), :), H, Gam, noise * eye(3));
    risks(c) = certaincoast.network.thresholdGoalEig(Sigma, -U(g.N + supply), combs(c, :), noise, 13, -limit);
end
[best, ib] = max(gains);
[bestRisk, ir] = max(risks);
fprintf('Salt mass near the toe: greedy design %s gives %.3f nats; exhaustive optimum %s %.3f; median design %.3f\n', ...
        mat2str(sort(S)), gain, mat2str(combs(ib, :)), best, median(gains));
fprintf('Supply cell above its limit (prior probability %.2f): best design %s gives %.3f nats;\n', ...
        0.5 * erfc(0.3 / sqrt(2)), mat2str(combs(ir, :)), bestRisk);
fprintf('  the salt-mass design %s gives it %.3f nats (%.0f%% of that optimum).\n', mat2str(combs(ib, :)), risks(ib), ...
        100 * risks(ib) / bestRisk);
