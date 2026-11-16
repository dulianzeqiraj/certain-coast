function F = fitMotionSeries(t, d, opts)
%FITMOTIONSERIES Rate, steps, seasonal cycle and post-seismic motion of a displacement series, with honest errors.
%   F = certaincoast.motion.fitMotionSeries(T, D) fits the displacement D [mm] observed at the times T [decimal years]
%   of one InSAR point or GNSS station by least squares with the model
%       d(t) = a + v (t - t0) + sum_k s_k H(t - T_k) + c1 sin(2 pi t) + c2 cos(2 pi t) + c3 sin(4 pi t) + c4 cos(4 pi t)
%              + p log(1 + (t - T_e)/tau) H(t - T_e),
%   with v the rate [mm/yr], s_k the steps at the epochs T_k (an earthquake, a change of equipment), the annual and
%   semi-annual cycle (the seasonal swelling and shrinking of soft ground with the water table) and an optional
%   logarithmic post-seismic term. F = certaincoast.motion.fitMotionSeries(T, D, OPTS) sets
%       steps        epochs of the steps [decimal years]                 ([])
%       seasonal     fit the annual and semi-annual cycle                (true)
%       postseismic  struct with fields epoch and tau [years], or []     ([])
%       noise        'ar1' (default): the errors of a displacement series are correlated in time, so the
%                    covariance of white noise is inflated by (1 + phi)/(1 - phi), phi the lag-one autocorrelation
%                    of the residuals (at most 0.95); 'white': no inflation
%       t0           reference epoch of the offset a                      (the mean of T)
%   Rates from white-noise least squares on InSAR or GNSS series are typically too confident by a factor of two or
%   more; the AR(1) inflation is a simple, conservative correction.
%
%   F has the fields rate and rateSigma [mm/yr], step and stepSigma [mm] (one value per epoch), annualAmplitude [mm],
%   post and postSigma [mm per unit log], offset, phi (lag-one autocorrelation), inflation, sigma0 (rms residual),
%   x and C (all coefficients and their covariance), names (the coefficient names), model and residual.
%
%   Example
%       t = 2016 + (0:6:3650).' / 365.25;  d = -4 * (t - 2016) + 2 * sin(2 * pi * t) + randn(size(t));
%       F = certaincoast.motion.fitMotionSeries(t, d);          % F.rate close to -4 mm/yr
%
%   See also certaincoast.motion.losToVertical, certaincoast.motion.relativeSeaLevel.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 3 || isempty(opts)
    opts = struct();
end
o = struct('steps', [], 'seasonal', true, 'postseismic', [], 'noise', 'ar1', 't0', []);
names = fieldnames(opts);
for i = 1:numel(names)
    if ~isfield(o, names{i})
        error('certaincoast:motion:option', 'unknown option %s.', names{i});
    end
    o.(names{i}) = opts.(names{i});
end
t = t(:);
d = d(:);
ok = isfinite(t) & isfinite(d);
t = t(ok);
d = d(ok);
if isempty(o.t0)
    o.t0 = mean(t);
end
n = numel(t);
A = [ones(n, 1), t - o.t0];
cols = {'offset', 'rate'};
steps = o.steps(:).';
for k = 1:numel(steps)
    A = [A, double(t >= steps(k))]; %#ok<AGROW>
    cols{end + 1} = sprintf('step %d', k); %#ok<AGROW>
end
if o.seasonal
    A = [A, sin(2 * pi * t), cos(2 * pi * t), sin(4 * pi * t), cos(4 * pi * t)];
    cols = [cols, {'annual sin', 'annual cos', 'semiannual sin', 'semiannual cos'}];
end
if ~isempty(o.postseismic)
    after = t >= o.postseismic.epoch;
    g = zeros(n, 1);
    g(after) = log(1 + (t(after) - o.postseismic.epoch) / o.postseismic.tau);
    A = [A, g];
    cols{end + 1} = 'postseismic';
end
p = size(A, 2);
if n <= p
    error('certaincoast:motion:tooShort', 'the series has %d values for %d coefficients.', n, p);
end
if rank(A) < p
    error('certaincoast:motion:rank', 'the model cannot be fitted: is a step or post-seismic epoch outside the series?');
end
x = A \ d;
r = d - A * x;
s2 = (r.' * r) / (n - p);
C = s2 * inv(A.' * A); %#ok<MINV>
phi = 0;
inflation = 1;
if strcmpi(o.noise, 'ar1')
    phi = min(max(sum(r(2:end) .* r(1:end - 1)) / sum(r .^ 2), 0), 0.95);
    inflation = (1 + phi) / (1 - phi);
elseif ~strcmpi(o.noise, 'white')
    error('certaincoast:motion:noise', 'noise must be ''ar1'' or ''white''.');
end
C = inflation * C;
C = 0.5 * (C + C.');
sx = sqrt(diag(C));
F = struct('rate', x(2), 'rateSigma', sx(2), 'offset', x(1), 'phi', phi, 'inflation', inflation, ...
           'sigma0', sqrt(s2));
F.step = x(3:2 + numel(steps)).';
F.stepSigma = sx(3:2 + numel(steps)).';
F.annualAmplitude = NaN;
if o.seasonal
    j = 3 + numel(steps);
    F.annualAmplitude = hypot(x(j), x(j + 1));
end
F.post = NaN;
F.postSigma = NaN;
if ~isempty(o.postseismic)
    F.post = x(end);
    F.postSigma = sx(end);
end
F.x = x;
F.C = C;
F.names = cols;
F.t = t;
F.model = A * x;
F.residual = r;
end
