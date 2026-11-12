function [risk, info] = rechargeRisk(A, cells, C, mu, hc)
%RECHARGERISK Where managed recharge lowers most the probability that a supply well is salinised (risk-oriented MAR).
%   [RISK, INFO] = certaincoast.network.rechargeRisk(A, CELLS, C, MU, HC) for the steady flow operator A (A h = rhs;
%   a recharge Q at cell k adds -Q e_k to rhs), one or two supply wells at CELLS, the covariance C and mean MU of
%   their heads, and the critical head HC (Ghyben-Herzberg: a well is intruded when h < HC) returns, for every cell,
%   the reduction of the intrusion probability per unit of recharge there (larger is better):
%
%       one well   P = Phi(-a),              -dP/dQ_k = phi(a) / s g_k
%       two wells  P = 1 - Phi2(a, b; r),    -dP/dQ_k = phi(a)/s_A Phi((b - r a)/sqrt(1 - r^2)) g_A,k + (A <-> B)
%
%   with a = (MU_A - HC)/s_A, b = (MU_B - HC)/s_B and g_w = -A^-T e_w (one transposed solve per well). A well that is
%   surely safe or surely lost gets no weight. INFO has fields P0 (prior intrusion probability), z, s, r, weights, g.
%
%   See also certaincoast.util.bvnCdf.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

cells = cells(:);
nw = numel(cells);
mu = mu(:);
if ~(nw == 1 || nw == 2) || ~isequal(size(C), [nw nw]) || numel(mu) ~= nw
    error('certaincoast:network:wells', 'One or two supply wells, with C nw-by-nw and mu nw-by-1.');
end
Nc = size(A, 1);
E = full(sparse(cells, (1:nw).', 1, Nc, nw));
g = -(A' \ E);
s = sqrt(diag(C));
z = (mu - hc) ./ s;
nphi = @(x) exp(-0.5 * x .^ 2) / sqrt(2 * pi);
ncdf = @(x) 0.5 * erfc(-x / sqrt(2));
if nw == 1
    r = NaN;
    w = nphi(z) / s;
    P0 = ncdf(-z);
else
    r = C(1, 2) / (s(1) * s(2));
    q = sqrt(1 - r ^ 2);
    w = [nphi(z(1)) * ncdf((z(2) - r * z(1)) / q) / s(1); nphi(z(2)) * ncdf((z(1) - r * z(2)) / q) / s(2)];
    P0 = 1 - certaincoast.util.bvnCdf(z(1), z(2), r);
end
risk = full(g * w);
info = struct('P0', P0, 'z', z.', 's', s.', 'r', r, 'weights', w.', 'g', full(g));
end
