function A = conductanceMatrix(n, p, sea)
%CONDUCTANCEMATRIX Finite-volume conductance matrix of a chain of N cells.
%   A = certaincoast.boundary.conductanceMatrix(N, P) is the N-by-N matrix of the conductances P.c between
%   neighbouring cells, with the first cell connected to the sea (conductance P.cs, head 0).
%   A = certaincoast.boundary.conductanceMatrix(N, P, false) leaves out the connection to the sea.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 3
    sea = true;
end
e = ones(n, 1);
A = p.c * (diag([e(1:end - 1); 0] + [0; e(2:end)]) - diag(e(1:end - 1), 1) - diag(e(1:end - 1), -1));
if sea
    A(1, 1) = A(1, 1) + p.cs;
end
end
