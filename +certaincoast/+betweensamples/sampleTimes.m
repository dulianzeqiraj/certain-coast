function t = sampleTimes(days, step, jitter)
%SAMPLETIMES Irregular sampling dates: a nominal step with a random jitter.
%   T = certaincoast.betweensamples.sampleTimes(DAYS, STEP) returns the sampling times (days, column vector) in
%   [0, DAYS) for a nominal step STEP with a jitter of +-20%.
%   T = certaincoast.betweensamples.sampleTimes(DAYS, STEP, JITTER) sets the relative jitter.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 3
    jitter = 0.2;
end
t = zeros(ceil(days / (step * (1 - jitter))) + 1, 1);
n = 0;
s = 0;
while s < days
    n = n + 1;
    t(n) = s;
    s = s + step * (1 + jitter * (2 * rand - 1));
end
t = t(1:n);
end
