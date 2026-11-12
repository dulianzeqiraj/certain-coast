%EXAMPLE_ALARM Calibrate a chloride alarm for a 10% miss rate and check it on new data.
%   Synthetic data only: no real well, asset or location is used. Run from the repository root.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

addpath(fileparts(fileparts(mfilename('fullpath'))));
rng(1);
threshold = 250;                                          % mg/L, drinking-water limit for chloride
n = 1200;
truth = 180 + 60 * sin(linspace(0, 12 * pi, n)).' + 15 * randn(n, 1);
forecast = truth + 20 * randn(n, 1);                      % an imperfect forecast

cal = 1:600;
live = 601:n;
margin = certaincoast.alarms.conformalMargin(truth(cal) - forecast(cal), 0.10);
alarm = certaincoast.alarms.issueAlarm(forecast(live), threshold, margin);
events = truth(live) >= threshold;
misses = sum(events & ~alarm);
bound = certaincoast.alarms.missRateUpperBound(misses, sum(events));
fprintf('margin %.1f mg/L; events %d, missed %d; 95%% upper bound on the miss rate %.3f\n', ...
        margin, sum(events), misses, bound);
