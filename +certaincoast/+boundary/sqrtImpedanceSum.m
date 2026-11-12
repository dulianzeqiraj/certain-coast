function [w, mu, M] = sqrtImpedanceSum(omegaMin, omegaMax, tol)
%SQRTIMPEDANCESUM Positive exponential sum for the half-order impedance of a long (semi-infinite) exterior.
%   [W, MU, M] = certaincoast.boundary.sqrtImpedanceSum(OMEGAMIN, OMEGAMAX, TOL) returns M compartments (weights W,
%   rates MU, all positive) whose impedance sum_m W(m) s / (s + MU(m)) approximates sqrt(s) for s = 1i*OMEGA,
%   OMEGA in [OMEGAMIN, OMEGAMAX], to a relative error of the order of TOL (Proposition 3 of the preparatory
%   experiment BEYOND-GO). The number of compartments grows with the logarithm of the frequency window and of 1/TOL:
%
%       M = ceil((X+ - X-) / ETA) + 1,   ETA = pi^2 / log(1/TOL),
%       X- = log(OMEGAMIN) + 2 log(pi TOL / 2),   X+ = log(OMEGAMAX) + 2 log(2 / (pi TOL)).
%
%   For a semi-infinite exterior of transmissivity T and storage S the impedance is sqrt(T S s); multiply W by
%   sqrt(T S).

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

L = log(1 / tol);
xm = log(omegaMin) + 2 * log(pi * tol / 2);
xp = log(omegaMax) + 2 * log(2 / (pi * tol));
M = ceil((xp - xm) / (pi ^ 2 / L)) + 1;
x = linspace(xm, xp, M).';
h = x(2) - x(1);
mu = exp(x);
w = (h / pi) * exp(x / 2);
end
