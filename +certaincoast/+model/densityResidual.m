function F = densityResidual(U, theta, g, P)
%DENSITYRESIDUAL Residual of the steady variable-density flow and salt transport in a vertical section.
%   F = certaincoast.model.densityResidual(U, THETA, G, P) for the state U = [h; c] (equivalent freshwater head and
%   relative salinity, 2 * G.N values), the log-conductivity THETA (G.N values, K = P.Kref * exp(THETA)), the grid G
%   (henryGrid) and the parameters P (henryParameters). In equivalent-freshwater-head form:
%
%       flow       div q = 0,                    q = -K (grad h + eps c e_z)
%       transport  div(q c) - div(phi Dmol grad c) = 0
%
%   Cell-centred finite volumes, harmonic-mean face conductivity, first-order upwind advection. Boundaries: inland
%   (x = 0) freshwater inflow P.Q with c = 0; sea (x = Lx) hydrostatic equivalent head with c = 1; base and top closed.
%   F has 2 * G.N values, flow first.
%
%   See also certaincoast.model.densityJacobian, certaincoast.model.solveDensity.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

Nx = g.Nx;
Ny = g.Ny;
N = g.N;
h = reshape(U(1:N), Nx, Ny).';                   % h(j, i): row j along y, column i along x
c = reshape(U(N + 1:2 * N), Nx, Ny).';
K = reshape(P.Kref * exp(theta(:)), Nx, Ny).';
hm = @(a, b) 2 * a .* b ./ (a + b);              % harmonic mean

qx = zeros(Ny, Nx + 1);
qy = zeros(Ny + 1, Nx);
qx(:, 2:Nx) = -hm(K(:, 1:Nx - 1), K(:, 2:Nx)) .* diff(h, 1, 2) / g.dx;
qy(2:Ny, :) = -hm(K(1:Ny - 1, :), K(2:Ny, :)) .* (diff(h, 1, 1) / g.dy + P.eps * 0.5 * (c(1:Ny - 1, :) + c(2:Ny, :)));
qx(:, 1) = P.Q;                                  % freshwater inflow
hSea = P.eps * (g.Ly - g.yc(:));                 % hydrostatic equivalent freshwater head of the sea
qx(:, Nx + 1) = -K(:, Nx) .* (hSea - h(:, Nx)) / (0.5 * g.dx);

Fh = diff(qx, 1, 2) / g.dx + diff(qy, 1, 1) / g.dy;

cW = [zeros(Ny, 1), c(:, 1:Nx - 1)];             % ghost values: fresh inland, sea water at the sea
cE = [c(:, 2:Nx), ones(Ny, 1)];
cS = [c(1, :); c(1:Ny - 1, :)];                  % closed base and top: zero gradient
cN = [c(2:Ny, :); c(Ny, :)];
qw = qx(:, 1:Nx);
qe = qx(:, 2:Nx + 1);
qs = qy(1:Ny, :);
qn = qy(2:Ny + 1, :);
up = @(cond, a, b) cond .* a + (~cond) .* b;     % upwind value
adv = (qe .* up(qe > 0, c, cE) - qw .* up(qw > 0, cW, c)) / g.dx ...
    + (qn .* up(qn > 0, c, cN) - qs .* up(qs > 0, cS, c)) / g.dy;
lap = (cE - 2 * c + cW) / g.dx ^ 2 + (cN - 2 * c + cS) / g.dy ^ 2;
Fc = adv - P.phi * P.Dmol * lap;
F = [reshape(Fh.', [], 1); reshape(Fc.', [], 1)];
end
