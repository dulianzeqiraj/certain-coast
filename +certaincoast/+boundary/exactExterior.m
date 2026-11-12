function [C, w, mu, Craw] = exactExterior(p)
%EXACTEXTERIOR The aquifer beyond the model boundary, replaced exactly by a Robin term and storage compartments.
%   [C, W, MU] = certaincoast.boundary.exactExterior(P) returns the representation of the exterior cells
%   P.nModel+1 ... P.nFull, seen from the last model cell, as a boundary flux
%
%       q = C u + sum_m W(m) (u - y_m),    dy_m/dt = MU(m) (u - y_m),
%
%   where u is the head of the boundary cell above its reference. All of C, W and MU are non-negative (Proposition 1
%   of the preparatory experiment BEYOND-GO): the boundary is passive, so it can never create energy, whatever the
%   forcing. The fourth output is C before it is clipped at zero (round-off for a closed exterior).
%
%   See also certaincoast.boundary.impedance, certaincoast.boundary.fitPositive.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

ne = p.nFull - p.nModel;
Ab = certaincoast.boundary.conductanceMatrix(ne, p, false);
Ab(1, 1) = Ab(1, 1) + p.c;                                  % coupling to the boundary cell
B = Ab / (p.S * p.dx);
[V, D] = eig((B + B') / 2);
mu = diag(D);
w = p.c ^ 2 * V(1, :).' .^ 2 ./ (p.S * p.dx * mu);
Craw = p.c - sum(w);
C = max(Craw, 0);
end
