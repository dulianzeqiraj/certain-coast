function p = aquifer(varargin)
%AQUIFER Parameters of the one-dimensional coastal aquifer used by the boundary functions.
%   P = certaincoast.boundary.aquifer() returns the default test aquifer: 10 km from the sea to a closed far end,
%   cells of 100 m, transmissivity 500 m2/d, storage coefficient 0.05, daily steps; the model sees only the first
%   2 km, and its inland boundary must stand in for the 8 km beyond.
%   P = certaincoast.boundary.aquifer('Name', Value, ...) overrides any field.
%
%   Fields
%       T, S, dx, dt     transmissivity (m2/d), storage coefficient, cell size (m), time step (d)
%       nFull, nModel    cells of the whole aquifer and of the model domain
%       well             cell of the pumping well; piezometer: an observation cell inside the model
%       bandCells        two observation cells next to the boundary (the last two before the boundary cell)
%       R0, Q0           mean recharge (m/d) and pumping (m2/d per metre of coast)
%       c, cs            conductance between cells and from the first cell to the sea (computed)

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

p = struct('T', 500, 'S', 0.05, 'dx', 100, 'dt', 1, 'nFull', 100, 'nModel', 20, 'well', 10, 'piezometer', 15, ...
           'bandCells', [18 19], 'R0', 5e-4, 'Q0', 0.6);
for i = 1:2:numel(varargin)
    p.(varargin{i}) = varargin{i + 1};
end
p.c = p.T / p.dx;
p.cs = 2 * p.T / p.dx;
p.boundaryCell = p.nModel;
end
