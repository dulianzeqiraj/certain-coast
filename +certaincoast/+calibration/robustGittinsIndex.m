function g = robustGittinsIndex(T, p, n)
%ROBUSTGITTINSINDEX Robust Gittins index at the state (p, n), read from a table of robustGittinsTable.
%   G = certaincoast.calibration.robustGittinsIndex(T, P, N) interpolates the table T bilinearly at the miss rates P
%   and the counts N (arrays of the same size, or one of them scalar). P is clipped to [0, 1] and N to the counts of
%   the table. After n plays with s misses an arm is read at P = (s + 1/2)/(n + 1) and N = n + 1.
%
%   Example
%       T = certaincoast.calibration.robustGittinsTable(0.5, 0.95, 60, 101, 20, 201);
%       g = certaincoast.calibration.robustGittinsIndex(T, [0.1 0.3], 20);
%
%   See also certaincoast.calibration.robustGittinsTable.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if ~isstruct(T) || ~all(isfield(T, {'p', 'n', 'gam'}))
    error('certaincoast:calibration:badTable', 'T must come from robustGittinsTable.');
end
if isscalar(p) && ~isscalar(n), p = p + zeros(size(n)); end
if isscalar(n) && ~isscalar(p), n = n + zeros(size(p)); end
p = min(max(p, 0), 1);
n = min(max(n, T.n(1)), T.n(end));
if numel(T.n) == 1
    g = reshape(interp1(T.p, T.gam(:, 1), p(:), 'linear'), size(p));
else
    g = interp2(T.n, T.p, T.gam, n, p, 'linear');
end
end
