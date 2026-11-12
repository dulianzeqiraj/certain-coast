function [C, w, mu] = fitPositive(Ztarget, nTerms, band, nFreq)
%FITPOSITIVE A reduced passive boundary: few compartments with non-negative weights fitted to a target impedance.
%   [C, W, MU] = certaincoast.boundary.fitPositive(ZTARGET, NTERMS, BAND) fits C >= 0 and W >= 0 for NTERMS rates MU,
%   spread geometrically from BAND(1)/3 to 3*BAND(2) (rad/day), so that impedance(1i*OMEGA, C, W, MU) matches
%   ZTARGET(1i*OMEGA) in relative terms over OMEGA in BAND. ZTARGET is a function handle of a column of complex
%   frequencies. The fit is a non-negative least-squares problem (lsqnonneg), so the reduced boundary is passive.
%   [...] = certaincoast.boundary.fitPositive(..., NFREQ) sets the number of fitting frequencies (default 240).
%
%   See also certaincoast.boundary.exactExterior, certaincoast.boundary.impedance.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 4
    nFreq = 240;
end
om = logspace(log10(band(1)), log10(band(2)), nFreq).';
mu = logspace(log10(band(1) / 3), log10(band(2) * 3), nTerms).';
s = 1i * om;
Zt = Ztarget(s);
A = [ones(size(s)), bsxfun(@rdivide, s, bsxfun(@plus, s, mu.'))];
A = bsxfun(@rdivide, A, abs(Zt));
b = Zt ./ abs(Zt);
x = lsqnonneg([real(A); imag(A)], [real(b); imag(b)]);
C = x(1);
w = x(2:end);
end
