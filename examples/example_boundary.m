%EXAMPLE_BOUNDARY The inland boundary of a coastal model: exact, reduced, plain, and learned from two piezometers.
%   A 10 km aquifer is modelled only over its first 2 km. The boundary at 2 km must stand in for the 8 km beyond.
%   The example compares, at the boundary cell, the exact passive boundary, a reduced one with 12 compartments, the
%   usual no-flow and fixed-head boundaries, and a boundary learned from two noisy piezometers, which then predicts
%   a drought with doubled pumping that the calibration never saw. Synthetic data only. Run from the repository root
%   (about 10 s).

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

addpath(fileparts(fileparts(mfilename('fullpath'))));
rng(7);
p = certaincoast.boundary.aquifer();
A = certaincoast.boundary.conductanceMatrix(p.nFull, p, true);
f = zeros(p.nFull, 1);
f(1:p.nModel) = p.R0 * p.dx;
f(p.well) = f(p.well) - p.Q0;
hs = A \ f;
H = hs(p.boundaryCell);
h0 = hs(1:p.nModel);
bc = p.boundaryCell;

nd = round(20 * 365.25);
[~, r, Q] = certaincoast.boundary.forcing(p, nd, 'broadband');
full = certaincoast.boundary.simulateFull(p, r, Q);
rmsBc = @(h, ref) sqrt(mean((h(:, bc) - ref(:, bc)) .^ 2));

[C, w, mu] = certaincoast.boundary.exactExterior(p);
[Cf, wf, muf] = certaincoast.boundary.fitPositive(@(s) certaincoast.boundary.impedance(s, C, w, mu), 12, ...
                                                 [2 * pi / (40 * 365.25), 2 * pi / 5]);
fprintf('20 years, error at the boundary cell (m, rms):\n');
fprintf('  exact passive boundary (%d compartments) %.1e\n', numel(mu), ...
        rmsBc(certaincoast.boundary.simulateTruncated(p, r, Q, C, w, mu, H, h0), full));
fprintf('  reduced passive boundary (12 compartments) %.4f\n', ...
        rmsBc(certaincoast.boundary.simulateTruncated(p, r, Q, Cf, wf, muf, H, h0), full));
fprintf('  no flow %.4f, fixed head %.4f\n', ...
        rmsBc(certaincoast.boundary.simulateTruncated(p, r, Q, 0, 0, 1, H, h0), full), ...
        rmsBc(certaincoast.boundary.simulateTruncated(p, r, Q, p.c, 0, 1, H, h0), full));

o1 = full(:, p.bandCells(1)) + 0.002 * randn(nd, 1);       % two piezometers with 2 mm noise
o2 = full(:, p.bandCells(2)) + 0.002 * randn(nd, 1);
[Ci, wi, mui] = certaincoast.boundary.identifyFromPiezometers(p, o1, o2, r, H);
[~, rp, Qp] = certaincoast.boundary.forcing(p, round(10 * 365.25), 'seasonal', [3 * 365.25, 8 * 365.25, 0.3], ...
                                            [2 * 365.25, Inf, 1.2]);
ref = certaincoast.boundary.simulateFull(p, rp, Qp);
fprintf('Prediction of a drought with doubled pumping, 10 years (m, rms at the boundary cell):\n');
fprintf('  boundary learned from two piezometers %.4f\n', ...
        rmsBc(certaincoast.boundary.simulateTruncated(p, rp, Qp, Ci, wi, mui, H, h0), ref));
fprintf('  no flow %.4f, fixed head %.4f\n', ...
        rmsBc(certaincoast.boundary.simulateTruncated(p, rp, Qp, 0, 0, 1, H, h0), ref), ...
        rmsBc(certaincoast.boundary.simulateTruncated(p, rp, Qp, p.c, 0, 1, H, h0), ref));
