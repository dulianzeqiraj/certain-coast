function [h, flux, energy] = simulateTruncated(p, r, Q, C, w, mu, H, h0, eta)
%SIMULATETRUNCATED The model domain only, closed by a passive learned boundary.
%   [H, FLUX, ENERGY] = certaincoast.boundary.simulateTruncated(P, R, Q, C, W, MU, HREF, H0) runs cells 1 ... P.nModel
%   by implicit Euler. The flux out of the boundary cell is
%
%       q = C u + sum_m W(m) (u - y_m),   dy_m/dt = MU(m) (u - y_m),   u = h(P.nModel) - HREF,
%
%   with the compartments y at rest at the start and the heads starting from H0 (column, P.nModel values).
%   FLUX is the boundary flux at each step, and ENERGY the stored energy of the heads (relative to HREF) and of the
%   compartments; with no forcing it never increases when C, W >= 0 (passivity).
%   [...] = certaincoast.boundary.simulateTruncated(..., ETA) adds a sea level ETA (m above its initial value).
%
%   Plain boundaries are special cases: no flow (C = 0, W = 0), fixed head (C = P.c, W = 0) and a general-head
%   boundary of conductance C (W = 0).

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

n = p.nModel;
bc = p.boundaryCell;
w = w(:);
mu = mu(:);
A = certaincoast.boundary.conductanceMatrix(n, p, true);
M = p.S * p.dx;
a = 1 ./ (1 + p.dt * mu);
G = C + sum(w .* a);                                       % implicit-Euler effective conductance
A(bc, bc) = A(bc, bc) + G;
[L, U, P] = lu(M / p.dt * eye(n) + A);
y = zeros(size(mu));
hk = h0(:);
h = zeros(numel(r), n);
flux = zeros(numel(r), 1);
energy = zeros(numel(r), 1);
store = zeros(size(mu));
store(mu > 0) = w(mu > 0) ./ mu(mu > 0);
for k = 1:numel(r)
    rhs = M / p.dt * hk;
    rhs(1:n) = rhs(1:n) + r(k) * p.dx;
    rhs(p.well) = rhs(p.well) - Q(k);
    rhs(bc) = rhs(bc) + G * H + sum(w .* a .* y);
    if nargin >= 9 && ~isempty(eta)
        rhs(1) = rhs(1) + p.cs * eta(k);
    end
    hk = U \ (L \ (P * rhs));
    u = hk(bc) - H;
    y = a .* (y + p.dt * mu * u);
    h(k, :) = hk.';
    flux(k) = C * u + sum(w .* (u - y));
    energy(k) = 0.5 * M * sum((hk - H) .^ 2) + 0.5 * sum(store .* y .^ 2);
end
end
