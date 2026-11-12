function m = conformalMargin(residuals, alpha)
%CONFORMALMARGIN One-sided split-conformal margin of an alarm.
%   M = certaincoast.alarms.conformalMargin(RESIDUALS, ALPHA) returns the margin M such that, for a new forecast error
%   exchangeable with the calibration errors, OBSERVED <= FORECAST + M holds with probability at least 1 - ALPHA.
%
%   RESIDUALS  calibration errors OBSERVED - FORECAST (positive when the forecast was too low), any shape.
%   ALPHA      target miss rate, in (0, 1).
%
%   M is the k-th smallest residual, k = ceil((n + 1)(1 - ALPHA)). It is Inf when there are too few calibration
%   errors for that level (k > n).
%
%   Example
%       m = certaincoast.alarms.conformalMargin(pastObserved - pastForecast, 0.10);   % 10% miss rate
%
%   See also certaincoast.alarms.issueAlarm, certaincoast.alarms.missRateUpperBound.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if ~isscalar(alpha) || ~isreal(alpha) || ~(alpha > 0 && alpha < 1)
    error('certaincoast:alarms:badAlpha', 'alpha must be a scalar in (0, 1).');
end
r = sort(double(residuals(:)));
if isempty(r) || any(~isfinite(r)) || ~isreal(r)
    error('certaincoast:alarms:badResiduals', 'residuals must be a non-empty array of finite real numbers.');
end
n = numel(r);
k = ceil((n + 1) * (1 - alpha));
if k <= n
    m = r(k);
else
    m = Inf;
end
end
