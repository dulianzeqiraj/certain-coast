function M = planAquifer(opts)
%PLANAQUIFER Synthetic plan-view confined coastal aquifer, for calibration and design tests (no field data).
%   M = certaincoast.model.planAquifer() builds a confined aquifer of 10 km inland by 20 km along the coast on
%   square cells of 500 m. The sea is the western edge: the cells of the first column hold the sea level (head 0);
%   the other edges are closed. The conductivity follows three synthetic facies (gravel channel belts, sandy gravel,
%   silty sand) that meander inland; the aquifer is divided into pilot zones and holds a coastal strip, two supply
%   wells A and B and a set of observation wells.
%   M = certaincoast.model.planAquifer(OPTS) changes the defaults given by the fields of OPTS:
%       Lx, Ly     size inland and along the coast [km]                          (10, 20)
%       dx         cell size [km]                                                 (0.5)
%       b          thickness [m]                                                  (30)
%       K0         conductivity of the three facies [m/d]                         ([300 150 50])
%       recharge   mean recharge [m/d]                                            (5e-4)
%       zones      pilot zones inland and along the coast                        ([3 4])
%       zoneSd     prior standard deviation of a zone log-multiplier              (0.5)
%       stripKm    width of the coastal strip [km]                                (1.5)
%       supply     supply wells A and B, one row [x y] each [km]                 ([3.25 14.75; 3.25 5.25])
%       nObs       number of observation wells                                    (25)
%       hc         critical head at a supply well [m]                             (1)
%
%   The critical head is the sharp-interface (Ghyben-Herzberg) reading: with the base of the aquifer 40 m below sea
%   level and a density ratio of 0.025, the toe of the salt wedge reaches a well where the steady head falls below
%   40 x 0.025 = 1 m.
%
%   M has the fields
%       nx, ny, dx, x, y      grid (cell centres x inland and y along the coast, km); X, Y: centres of all Nc cells
%       Nc, b, recharge, cellArea (m^2), hc
%       K0, K0c, facies       facies conductivities, background conductivity of every cell [m/d], facies class 1-3
%       zone, Nz, zoneFacies  pilot zone of every cell, number of zones, main facies class of every zone
%       zonePrior             prior covariance of the zone log-multipliers: zoneSd^2 (I + E)/2, E(i,j) = 1 when
%                             zones i and j share a facies class (positive definite)
%       coastal, strip        cells at the sea (head 0) and within stripKm of the coast (logical, Nc x 1);
%                             stripKm, the width of the strip [km]
%       supplyCell, obsCell   cells of the supply wells A and B and of the observation wells
%       interior, pos, ek, el, faceInterior   cells off the coast and their position, the faces (cell pairs) and
%                                             the faces between two interior cells, used by planHeads
%   The observation wells are the points of a Halton sequence (bases 2 and 3) at least 0.5 km from the coast and
%   from the supply wells, one per cell, so the layout is the same in MATLAB and in GNU Octave.
%
%   Example
%       M = certaincoast.model.planAquifer();
%       H = certaincoast.model.planHeads(M, zeros(M.Nc, 1), M.recharge);
%       fprintf('natural head at A %.2f m, at B %.2f m\n', H(M.supplyCell, 1));
%
%   See also certaincoast.model.planHeads, certaincoast.calibration.trialGo.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 1 || isempty(opts)
    opts = struct();
end
o = struct('Lx', 10, 'Ly', 20, 'dx', 0.5, 'b', 30, 'K0', [300 150 50], 'recharge', 5e-4, 'zones', [3 4], ...
           'zoneSd', 0.5, 'stripKm', 1.5, 'supply', [3.25 14.75; 3.25 5.25], 'nObs', 25, 'hc', 1);
names = fieldnames(opts);
for i = 1:numel(names)
    if ~isfield(o, names{i})
        error('certaincoast:model:option', 'unknown option %s.', names{i});
    end
    o.(names{i}) = opts.(names{i});
end
if numel(o.K0) ~= 3 || any(o.K0(:) <= 0) || o.b <= 0 || o.dx <= 0 || o.recharge < 0
    error('certaincoast:model:badAquifer', 'need three positive facies conductivities, b > 0, dx > 0, recharge >= 0.');
end
if ~isequal(size(o.supply), [2 2])
    error('certaincoast:model:badSupply', 'supply must hold two rows [x y] (km).');
end

nx = round(o.Lx / o.dx);
ny = round(o.Ly / o.dx);
x = ((1:nx) - 0.5) * o.dx;
y = ((1:ny) - 0.5) * o.dx;
[XX, YY] = meshgrid(x, y);                         % ny x nx, cell (j, i) has the index j + (i - 1) ny
X = XX(:);
Y = YY(:);
Nc = nx * ny;
idx = reshape(1:Nc, ny, nx);

% three synthetic facies: gravel channel belts that meander inland, sandy gravel, silty sand
s = sin(2 * pi * (Y - 0.9 * sin(2 * pi * X / 5)) / 3.3) + 0.35 * cos(2 * pi * X / o.Lx);
facies = 2 * ones(Nc, 1);
facies(s > 0.5) = 1;
facies(s < -0.5) = 3;
K0c = reshape(o.K0(facies), [], 1);

% pilot zones (blocks), their main facies and the facies prior of their log-multipliers
zx = min(max(ceil(X / (o.Lx / o.zones(1))), 1), o.zones(1));
zy = min(max(ceil(Y / (o.Ly / o.zones(2))), 1), o.zones(2));
zone = (zx - 1) * o.zones(2) + zy;
Nz = prod(o.zones);
zoneFacies = zeros(Nz, 1);
for z = 1:Nz
    zoneFacies(z) = mode(facies(zone == z));
end
zonePrior = o.zoneSd ^ 2 * 0.5 * (eye(Nz) + double(zoneFacies == zoneFacies.'));

coastal = X < o.dx;                                % first column: the sea
strip = X <= o.stripKm + 1e-9;

% faces between neighbouring cells, each once
ek = [reshape(idx(:, 1:end - 1), [], 1); reshape(idx(1:end - 1, :), [], 1)];
el = [reshape(idx(:, 2:end), [], 1); reshape(idx(2:end, :), [], 1)];
interior = find(~coastal);
pos = zeros(Nc, 1);
pos(interior) = 1:numel(interior);
faceInterior = ~coastal(ek) & ~coastal(el);

% supply wells: the nearest cell centres
supplyCell = zeros(1, 2);
for w = 1:2
    [~, i] = min(abs(x - o.supply(w, 1)));
    [~, j] = min(abs(y - o.supply(w, 2)));
    supplyCell(w) = idx(j, i);
end
if any(coastal(supplyCell))
    error('certaincoast:model:badSupply', 'a supply well lies at the sea.');
end

% observation wells: Halton points, at least 0.5 km from the coast and from the supply wells, one per cell
obsCell = zeros(1, 0);
k = 0;
while numel(obsCell) < o.nObs && k < 100 * o.nObs
    k = k + 1;
    px = 0.5 + (o.Lx - 0.75) * radicalInverse(k, 2);
    py = 0.25 + (o.Ly - 0.5) * radicalInverse(k, 3);
    [~, i] = min(abs(x - px));
    [~, j] = min(abs(y - py));
    c = idx(j, i);
    far = min(hypot(X(supplyCell) - X(c), Y(supplyCell) - Y(c))) >= 0.5;
    if ~coastal(c) && X(c) >= 0.5 && far && ~any(obsCell == c)
        obsCell(end + 1) = c; %#ok<AGROW>
    end
end

M = struct('nx', nx, 'ny', ny, 'dx', o.dx, 'x', x, 'y', y, 'X', X, 'Y', Y, 'Nc', Nc, 'b', o.b, ...
           'recharge', o.recharge, 'cellArea', (1000 * o.dx) ^ 2, 'hc', o.hc, 'K0', o.K0(:).', 'K0c', K0c, ...
           'facies', facies, 'zone', zone, 'Nz', Nz, 'zoneFacies', zoneFacies);
M.zonePrior = zonePrior;
M.coastal = coastal;
M.strip = strip;
M.stripKm = o.stripKm;
M.supplyCell = supplyCell;
M.obsCell = obsCell;
M.interior = interior;
M.pos = pos;
M.ek = ek;
M.el = el;
M.faceInterior = faceInterior;
end

function v = radicalInverse(i, base)
% van der Corput radical inverse of the integer i in the given base
v = 0;
f = 1 / base;
while i > 0
    v = v + f * mod(i, base);
    i = floor(i / base);
    f = f / base;
end
end
