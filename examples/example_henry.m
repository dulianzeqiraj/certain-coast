%EXAMPLE_HENRY The Henry seawater-intrusion benchmark with the variable-density model, on four grids.
%   Solves the steady coupled flow and salt transport of the Henry problem and prints the position of the toe (the
%   0.5 isochlor at the base) and the residual on grids from 16 x 8 to 61 x 31. Run from the repository root
%   (a few seconds in MATLAB, about half a minute in GNU Octave).

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

addpath(fileparts(fileparts(mfilename('fullpath'))));
P = certaincoast.model.henryParameters();
fprintf('Henry benchmark: Henry number phi*Dmol/Q = %.3f\n', P.phi * P.Dmol / P.Q);
for s = [16 8; 32 16; 41 21; 61 31].'
    g = certaincoast.model.henryGrid(s(1), s(2));
    t0 = tic;
    [U, normF, it] = certaincoast.model.solveDensity(zeros(g.N, 1), g, P, 1e-12);
    toe = certaincoast.model.toePosition(U, g);
    fprintf('  %2d x %2d: toe %.4f m (%.3f of the length), residual %.1e, %d Newton steps, %.1f s\n', ...
            s(1), s(2), toe, toe / g.Lx, normF, it, toc(t0));
end
