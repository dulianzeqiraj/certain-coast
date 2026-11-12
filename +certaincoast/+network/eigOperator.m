function [gain, IF, sigma0sq] = eigOperator(G, H, GamF, Gn)
%EIGOPERATOR Goal-oriented information operator and expected information gain of a monitoring design.
%   [GAIN, IF, SIGMA0SQ] = certaincoast.network.eigOperator(G, H, GAMF, GN) for observation sensitivities G (m-by-n,
%   one row per sensor of the design), goal sensitivities H (p-by-n; p = 1 for a scalar goal such as the salt mass
%   near the toe), prior covariance GAMF (n-by-n) of the uncertain parameters and noise covariance GN (m-by-m):
%
%       IF   = H GAMF G' (G GAMF G' + GN)^-1 G GAMF H'          (information operator)
%       GAIN = 1/2 log det(H GAMF H') / det(H GAMF H' - IF)       (nats; scalar goal: -1/2 log(1 - IF/SIGMA0SQ))
%
%   GAIN is the mutual information between the data and the goal (exact for a linear-Gaussian model, the Laplace
%   approximation otherwise), so sensors are valued by what they say about the decision, not about every parameter.
%
%   See also certaincoast.network.greedyDesign, certaincoast.network.adjointSensitivities.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

[m, n] = size(G);
[p, nH] = size(H);
if nH ~= n || ~isequal(size(GamF), [n n]) || ~isequal(size(Gn), [m m])
    error('certaincoast:network:dim', 'G is m-by-n, H p-by-n, GamF n-by-n and Gn m-by-m.');
end
GH = GamF * H';
sigma0sq = H * GH;
GgH = G * GH;
IF = GgH' * ((G * GamF * G' + Gn) \ GgH);
IF = 0.5 * (IF + IF');
if p == 1
    gain = -0.5 * log(max(1 - IF / sigma0sq, eps));
else
    gain = 0.5 * (logdetSpd(sigma0sq) - logdetSpd(sigma0sq - IF));
end
gain = max(gain, 0);
end

function ld = logdetSpd(A)
A = 0.5 * (A + A');
[R, flag] = chol(A);
if flag ~= 0
    A = A + (eps + abs(min(svd(A)))) * eye(size(A));
    R = chol(A);
end
ld = 2 * sum(log(diag(R)));
end
