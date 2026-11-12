function [J, F0, Jtheta] = densityJacobian(U, theta, g, P, epsFd)
%DENSITYJACOBIAN Sparse Jacobians of the variable-density residual by coloured finite differences.
%   [J, F0] = certaincoast.model.densityJacobian(U, THETA, G, P) returns the sparse Jacobian J = dF/dU and the
%   residual F0 at U. [J, F0, JTHETA] also returns dF/dTHETA, the derivative with respect to the log-conductivity.
%   [...] = certaincoast.model.densityJacobian(..., EPSFD) sets the relative difference step (default 1e-7).
%
%   Each residual depends only on the cell and its four neighbours, so cells whose colours mod(i + 5 j, 13) are equal
%   never share a residual: 13 colours per field give the whole Jacobian with 2 * 13 + 1 residual evaluations instead
%   of 2 * G.N + 1, and its derivative with respect to THETA with 13 more.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 5
    epsFd = 1e-7;
end
N = g.N;
[I, Jc] = ndgrid(1:g.Nx, 1:g.Ny);                % column i and row j of every cell, flat order k = (j-1) Nx + i
colour = mod(I(:) + 5 * Jc(:), 13);
F0 = certaincoast.model.densityResidual(U, theta, g, P);
nb = neighbours(g);
J = assemble(@(d, k) certaincoast.model.densityResidual(U + d, theta, g, P), U, 2 * N, [0, N], colour, nb, F0, epsFd);
if nargout > 2
    Jtheta = assemble(@(d, k) certaincoast.model.densityResidual(U, theta + d, g, P), theta(:), N, 0, colour, nb, ...
                      F0, epsFd);
end
end

function A = assemble(res, x, n, offsets, colour, nb, F0, epsFd)
N = numel(colour);
rows = [];
cols = [];
vals = [];
for off = offsets
    for k = 0:12
        cells = find(colour == k);
        idx = off + cells;
        step = epsFd * max(1, abs(x(idx)));
        d = zeros(n, 1);
        d(idx) = step;
        dF = res(d, k) - F0;
        for m = 1:numel(cells)
            r = [nb{cells(m)}; N + nb{cells(m)}];
            rows = [rows; r]; %#ok<AGROW>
            cols = [cols; idx(m) * ones(numel(r), 1)]; %#ok<AGROW>
            vals = [vals; dF(r) / step(m)]; %#ok<AGROW>
        end
    end
end
A = sparse(rows, cols, vals, 2 * N, n);
end

function nb = neighbours(g)
nb = cell(g.N, 1);
for j = 1:g.Ny
    for i = 1:g.Nx
        list = (j - 1) * g.Nx + i;
        if i > 1, list(end + 1) = list(1) - 1; end %#ok<AGROW>
        if i < g.Nx, list(end + 1) = list(1) + 1; end %#ok<AGROW>
        if j > 1, list(end + 1) = list(1) - g.Nx; end %#ok<AGROW>
        if j < g.Ny, list(end + 1) = list(1) + g.Nx; end %#ok<AGROW>
        nb{list(1)} = list(:);
    end
end
end
