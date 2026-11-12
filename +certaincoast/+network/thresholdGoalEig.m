function [gain, p0] = thresholdGoalEig(Sigma, mu, design, noiseVar, goal, threshold)
%THRESHOLDGOALEIG Information a design gives about a yes/no question: is a supply well salinised? (PRAG-GO)
%   [GAIN, P0] = certaincoast.network.thresholdGoalEig(SIGMA, MU, DESIGN, NOISEVAR, GOAL, THRESHOLD) for jointly
%   Gaussian variables with covariance SIGMA (candidate observations and goal heads together), the means MU of the
%   goal heads, the observed variables DESIGN (indices into SIGMA, read with noise variance NOISEVAR), one or two goal
%   heads GOAL (indices into SIGMA) and their THRESHOLD (scalar or one per goal), returns the mutual information (nats)
%   between the data and the binary event
%
%       one goal head    h_A < THRESHOLD                  (well A intruded)
%       two goal heads   h_A < THRESHOLD(1) or h_B < THRESHOLD(2)   (A or B intruded)
%
%   GAIN = H(P0) - E[H(P | data)], H the binary entropy, by Gauss-Hermite quadrature (24 nodes per dimension) over the
%   possible data. P0 is the prior probability of the event. Unlike the information about a continuous quantity, it
%   values a sensor by how much it can change the answer to the question the operator must decide.
%
%   See also certaincoast.network.exactGoalEig, certaincoast.network.eigOperator.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

goal = goal(:);
mu = mu(:);
t = threshold(:) .* ones(numel(goal), 1);
[hx, hw] = hermite(24);
Sy = Sigma(design, design) + noiseVar * eye(numel(design));
Cz = Sigma(goal, design);
K = Cz * (Sy \ Cz');                             % covariance explained by the data
Sz = Sigma(goal, goal);
if numel(goal) == 1
    R2 = min(K / Sz, 1 - 1e-12);
    a = (t - mu) / sqrt(Sz);
    p0 = Phi(a);
    p = Phi((a - sqrt(R2) * hx) / sqrt(1 - R2));
    gain = hbin(p0) - sum(hw .* hbin(p));
    return
end
Pst = Sz - K;                                    % posterior covariance of the two goal heads
K = 0.5 * (K + K');
[V, D] = eig(K);
L = V * diag(sqrt(max(diag(D), 0)));
[U1, U2] = ndgrid(hx, hx);
shift = L * [U1(:).'; U2(:).'];
W = reshape(hw * hw.', 1, []);
s = sqrt(max(diag(Pst), 1e-12));
r = min(max(Pst(1, 2) / (s(1) * s(2)), -0.999), 0.999);
p = 1 - certaincoast.util.bvnCdf((mu(1) + shift(1, :) - t(1)) / s(1), (mu(2) + shift(2, :) - t(2)) / s(2), r);
r0 = Sz(1, 2) / sqrt(Sz(1, 1) * Sz(2, 2));
p0 = 1 - certaincoast.util.bvnCdf((mu(1) - t(1)) / sqrt(Sz(1, 1)), (mu(2) - t(2)) / sqrt(Sz(2, 2)), r0);
gain = hbin(p0) - W * hbin(p(:));
end

function v = Phi(x)
v = 0.5 * erfc(-x / sqrt(2));
end

function v = hbin(p)
p = min(max(p, 0), 1);
v = zeros(size(p));
k = p > 0 & p < 1;
v(k) = -(p(k) .* log(p(k)) + (1 - p(k)) .* log(1 - p(k)));
end

function [x, w] = hermite(n)
% probabilists' Gauss-Hermite nodes and normalised weights (Golub-Welsch)
J = diag(sqrt(1:n - 1), 1) + diag(sqrt(1:n - 1), -1);
[V, D] = eig(J);
[x, o] = sort(diag(D));
w = V(1, o).' .^ 2;
w = w / sum(w);
end
