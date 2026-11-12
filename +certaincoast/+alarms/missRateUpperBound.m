function u = missRateUpperBound(misses, events, confidence)
%MISSRATEUPPERBOUND Exact (Clopper-Pearson) upper confidence bound on the miss rate observed in operation.
%   U = certaincoast.alarms.missRateUpperBound(MISSES, EVENTS) is the 95% upper bound on the probability of missing a
%   dangerous event after MISSES misses among EVENTS dangerous events.
%   U = certaincoast.alarms.missRateUpperBound(MISSES, EVENTS, CONFIDENCE) uses the given confidence level.
%
%   The bound is the CONFIDENCE quantile of the Beta(MISSES + 1, EVENTS - MISSES) law; it is 1 when EVENTS = 0 or
%   MISSES = EVENTS. Only core functions are used (betaincinv), so no toolbox is needed.
%
%   Example
%       u = certaincoast.alarms.missRateUpperBound(1, 20);   % 0.216
%
%   See also certaincoast.alarms.conformalMargin.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 3
    confidence = 0.95;
end
if ~isscalar(misses) || ~isscalar(events) || misses ~= round(misses) || events ~= round(events) ...
        || events < 0 || misses < 0 || misses > events
    error('certaincoast:alarms:badCounts', 'need integers 0 <= misses <= events.');
end
if ~isscalar(confidence) || ~(confidence > 0 && confidence < 1)
    error('certaincoast:alarms:badConfidence', 'confidence must be a scalar in (0, 1).');
end
if events == 0 || misses == events
    u = 1;
else
    u = betaincinv(confidence, misses + 1, events - misses);
end
end
