function h = simulateFull(p, r, Q, eta)
%SIMULATEFULL The whole aquifer, from the sea to its far end: the reference the boundary must reproduce.
%   H = certaincoast.boundary.simulateFull(P, R, Q) runs the full aquifer (P.nFull cells, closed far end) by implicit
%   Euler from its steady state for the mean recharge P.R0 and pumping P.Q0. R and Q are daily series; recharge acts
%   on the model domain (cells 1 ... P.nModel), pumping at cell P.well. H is numel(R)-by-P.nFull (m above sea level).
%   H = certaincoast.boundary.simulateFull(P, R, Q, ETA) adds a sea level ETA (m above its initial value).

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

n = p.nFull;
A = certaincoast.boundary.conductanceMatrix(n, p, true);
M = p.S * p.dx;
h = A \ sourceVector(n, p, p.R0, p.Q0);
[L, U, P] = lu(M / p.dt * eye(n) + A);
out = zeros(numel(r), n);
for k = 1:numel(r)
    f = sourceVector(n, p, r(k), Q(k));
    if nargin >= 4 && ~isempty(eta)
        f(1) = f(1) + p.cs * eta(k);
    end
    h = U \ (L \ (P * (M / p.dt * h + f)));
    out(k, :) = h.';
end
h = out;
end

function f = sourceVector(n, p, r, Q)
f = zeros(n, 1);
f(1:p.nModel) = r * p.dx;
f(p.well) = f(p.well) - Q;
end
