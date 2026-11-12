function g = henryGrid(Nx, Ny, Lx, Ly)
%HENRYGRID Cell-centred rectangular grid of a vertical coastal section.
%   G = certaincoast.model.henryGrid(NX, NY) is the grid of the Henry domain, 2 m long and 1 m high, with NX cells
%   along x (from the inland boundary, x = 0, to the sea, x = LX) and NY cells along y (upwards from the base).
%   G = certaincoast.model.henryGrid(NX, NY, LX, LY) sets the size of the section (m).
%   The state of cell (j, i) (row j along y, column i along x) has the flat index k = (j - 1) * NX + i.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 3
    Lx = 2.0;
end
if nargin < 4
    Ly = 1.0;
end
g = struct('Nx', Nx, 'Ny', Ny, 'Lx', Lx, 'Ly', Ly, 'dx', Lx / Nx, 'dy', Ly / Ny, 'N', Nx * Ny);
g.xc = ((1:Nx) - 0.5) * g.dx;
g.yc = ((1:Ny) - 0.5) * g.dy;
end
