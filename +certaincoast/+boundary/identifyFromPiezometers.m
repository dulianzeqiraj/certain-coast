function [C, w, mu, info] = identifyFromPiezometers(p, obs1, obs2, r, H, mu)
%IDENTIFYFROMPIEZOMETERS Learn the passive inland boundary from two piezometers next to it.
%   [C, W, MU] = certaincoast.boundary.identifyFromPiezometers(P, OBS1, OBS2, R, HREF) fits the boundary
%   q = C u + sum_m W(m) (u - y_m) (C >= 0, W >= 0) from daily heads OBS1 and OBS2 observed in the two cells
%   P.bandCells next to the boundary cell, and the daily recharge R. The head at the first piezometer drives a small
%   model of the last two cells, and the parameters are chosen so that the head it computes at the second piezometer
%   matches OBS2 (output error: the noisy heads are never differentiated). HREF is the reference head of the
%   boundary cell. The rates MU are fixed on a geometric grid (12 values from 5 days to 120 years by default).
%   [...] = certaincoast.boundary.identifyFromPiezometers(..., MU) uses the given rates.
%   [C, W, MU, INFO] also returns the fit diagnostics of certaincoast.util.boundedLeastSquares.
%
%   The two piezometers must lie in a band with no pumping (Proposition 4 of the preparatory experiment BEYOND-GO);
%   a varied forcing (pumping that changes over months and years) identifies the boundary over more time scales than
%   a purely seasonal one.
%
%   See also certaincoast.boundary.simulateTruncated, certaincoast.util.boundedLeastSquares.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 6 || isempty(mu)
    mu = logspace(log10(2 * pi / (40 * 365.25) / 3), log10(2 * pi / 5 * 3), 12).';
end
mu = mu(:);
obs1 = obs1(:);
obs2 = obs2(:);
r = r(:);
theta0 = p.c / (numel(mu) + 1) * ones(numel(mu) + 1, 1);
fun = @(theta) bandModel(theta, p, obs1, r, obs2(1), H, mu) - repmat(obs2, 1, size(theta, 2));
[theta, info] = certaincoast.util.boundedLeastSquares(fun, theta0, zeros(size(theta0)));
C = theta(1);
w = theta(2:end);
end

function out = bandModel(theta, p, o1, r, h1init, H, mu)
% Heads of the second piezometer for each column of THETA = [C; W] (vectorised over parameter sets).
M = p.S * p.dx;
c = p.c;
dt = p.dt;
a = 1 ./ (1 + dt * mu);
K = size(theta, 2);
w = theta(2:end, :);
G = theta(1, :) + sum(bsxfun(@times, w, a), 1);
wa = bsxfun(@times, w, a);
A11 = M / dt + 2 * c;
A12 = -c;
A21 = -c;
A22 = M / dt + c + G;
dt_ = A11 * A22 - A12 * A21;
h1 = h1init * ones(1, K);
h2 = H * ones(1, K);
y = zeros(numel(mu), K);
out = zeros(numel(o1), K);
for k = 1:numel(o1)
    b1 = M / dt * h1 + c * o1(k) + r(k) * p.dx;
    b2 = M / dt * h2 + G * H + sum(wa .* y, 1) + r(k) * p.dx;
    h1 = (b1 .* A22 - A12 * b2) ./ dt_;
    h2 = (A11 * b2 - A21 * b1) ./ dt_;
    y = bsxfun(@times, a, y + dt * bsxfun(@times, mu, h2 - H));
    out(k, :) = h1;
end
end
