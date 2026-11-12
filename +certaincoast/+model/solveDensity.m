function [U, normF, iterations] = solveDensity(theta, g, P, tol, maxIter)
%SOLVEDENSITY Steady variable-density flow and salt transport by damped Newton iterations.
%   [U, NORMF, ITERATIONS] = certaincoast.model.solveDensity(THETA, G, P) solves densityResidual(U) = 0 from a
%   uniform start (h = 0, c = 0.1), with the sparse coloured Jacobian and a backtracking line search, to a residual
%   norm of 1e-10. U = [h; c].
%   [...] = certaincoast.model.solveDensity(..., TOL, MAXITER) sets the tolerance and the iteration limit (60).
%
%   Example (Henry benchmark)
%       g = certaincoast.model.henryGrid(41, 21);
%       U = certaincoast.model.solveDensity(zeros(g.N, 1), g, certaincoast.model.henryParameters());
%       toe = certaincoast.model.toePosition(U, g)        % 1.392 m
%
%   See also certaincoast.model.densityResidual, certaincoast.model.toePosition.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 4 || isempty(tol)
    tol = 1e-10;
end
if nargin < 5 || isempty(maxIter)
    maxIter = 60;
end
N = g.N;
U = [zeros(N, 1); 0.1 * ones(N, 1)];
normF = Inf;
for iterations = 1:maxIter
    [J, F] = certaincoast.model.densityJacobian(U, theta, g, P);
    normF = norm(F);
    if normF < tol
        break
    end
    dU = -(J \ F);
    lambda = 1;
    for s = 1:20
        if norm(certaincoast.model.densityResidual(U + lambda * dU, theta, g, P)) < normF
            break
        end
        lambda = lambda / 2;
    end
    U = U + lambda * dU;
end
end
