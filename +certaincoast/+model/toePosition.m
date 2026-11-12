function toe = toePosition(U, g, level)
%TOEPOSITION Position of the saltwater toe: where the isochlor reaches the base of the aquifer.
%   TOE = certaincoast.model.toePosition(U, G) is the distance x (m, from the inland boundary) at which the relative
%   salinity of the bottom row of cells crosses 0.5, by linear interpolation between cell centres; NaN when it does not.
%   TOE = certaincoast.model.toePosition(U, G, LEVEL) uses another isochlor.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 3
    level = 0.5;
end
cb = U(g.N + 1:g.N + g.Nx);                      % bottom row (j = 1)
toe = NaN;
for i = 1:g.Nx - 1
    if (cb(i) - level) * (cb(i + 1) - level) < 0
        toe = g.xc(i) + (level - cb(i)) / (cb(i + 1) - cb(i)) * (g.xc(i + 1) - g.xc(i));
        return
    end
end
end
