function [vUp, vEast, sUp, sEast] = losToVertical(losA, losD, uA, uD, sA, sD)
%LOSTOVERTICAL Vertical and east-west motion from InSAR line-of-sight motion seen from two orbits.
%   [VUP, VEAST] = certaincoast.motion.losToVertical(LOSA, LOSD, UA, UD) combines the line-of-sight velocities (or
%   displacements) LOSA and LOSD of the same ground points seen from an ascending and a descending orbit, positive
%   towards the satellite, with the unit vectors UA and UD = [east north up] that point from the ground to each
%   satellite (as the EGMS products give them, for example). With the north motion neglected, to which the near-polar
%   orbits of Sentinel-1 are almost blind, each point solves
%       [UA(1) UA(3); UD(1) UD(3)] [VEAST; VUP] = [LOSA; LOSD].
%   UA and UD are 3-element vectors (one geometry for all points) or n-by-3 matrices (one row per point).
%   [VUP, VEAST, SUP, SEAST] = certaincoast.motion.losToVertical(..., SA, SD) also propagates the standard deviations
%   SA and SD of the two line-of-sight values, taken as independent.
%   Vertical land motion turns into relative sea level with certaincoast.motion.relativeSeaLevel.
%
%   Example (incidence about 39 degrees, ascending heading about -13 degrees, descending about -167 degrees)
%       uA = [-0.61 -0.14 0.78];  uD = [0.61 -0.14 0.78];
%       [vUp, vEast] = certaincoast.motion.losToVertical(-3.6, -2.1, uA, uD);
%
%   See also certaincoast.motion.fitMotionSeries, certaincoast.motion.relativeSeaLevel.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if ~isequal(size(losA), size(losD))
    error('certaincoast:motion:size', 'LOSA and LOSD must have the same size.');
end
n = numel(losA);
uA = geometry(uA, n);
uD = geometry(uD, n);
det = uA(:, 1) .* uD(:, 3) - uA(:, 3) .* uD(:, 1);
if any(abs(det) < 1e-3)
    error('certaincoast:motion:geometry', 'the two looks are too similar to separate east and vertical motion.');
end
a = losA(:);
b = losD(:);
vEast = (uD(:, 3) .* a - uA(:, 3) .* b) ./ det;
vUp = (uA(:, 1) .* b - uD(:, 1) .* a) ./ det;
vEast = reshape(vEast, size(losA));
vUp = reshape(vUp, size(losA));
if nargout > 2
    if nargin < 6
        error('certaincoast:motion:sigma', 'give SA and SD to get the standard deviations.');
    end
    sa = sA(:) + zeros(n, 1);
    sb = sD(:) + zeros(n, 1);
    sEast = sqrt((uD(:, 3) .* sa) .^ 2 + (uA(:, 3) .* sb) .^ 2) ./ abs(det);
    sUp = sqrt((uD(:, 1) .* sa) .^ 2 + (uA(:, 1) .* sb) .^ 2) ./ abs(det);
    sEast = reshape(sEast, size(losA));
    sUp = reshape(sUp, size(losA));
end
end

function u = geometry(u, n)
% one geometry for all points, or one row per point, each a unit vector [east north up]
if numel(u) == 3
    u = repmat(u(:).', n, 1);
end
if ~isequal(size(u), [n 3])
    error('certaincoast:motion:geometry', 'unit vectors must be 3-element vectors or n-by-3 matrices.');
end
if any(abs(sqrt(sum(u .^ 2, 2)) - 1) > 1e-2)
    error('certaincoast:motion:geometry', 'the line-of-sight vectors must have unit length.');
end
end
