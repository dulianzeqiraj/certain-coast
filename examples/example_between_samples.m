%EXAMPLE_BETWEEN_SAMPLES Does chloride cross 250 mg/L between two laboratory samples, and what must be measured?
%   Chloride with mean 180 mg/L and standard deviation 35 mg/L (limit 250 mg/L, so z = 2) and a correlation time of
%   8 days, over a season of 120 days. Records of 8 years with quarterly, monthly or weekly samples, or a daily
%   logger. For each design: the time scale the record allows, and the band of the probability of crossing the limit
%   in the season. Synthetic data only. Run from the repository root (about 10 s).

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

addpath(fileparts(fileparts(mfilename('fullpath'))));
rng(3);
z = (250 - 180) / 35;
kGrid = exp(linspace(log(1 / 400), log(4), 41));
curve = certaincoast.betweensamples.crossingProbability(kGrid, z, 120, 0.25, 8000);
truth = interp1(log(kGrid), curve, log(1 / 8));
fprintf('probability of crossing in the season %.2f; at the samples alone (marginal tail) %.3f\n', ...
        truth, 0.5 * erfc(z / sqrt(2)));
designs = {'quarterly samples', 91; 'monthly samples', 30; 'weekly samples', 7; 'daily logger', 1};
for d = 1:size(designs, 1)
    t = certaincoast.betweensamples.sampleTimes(8 * 365, designs{d, 2});
    x = certaincoast.betweensamples.ouAtTimes(t, 1 / 8);
    [kHat, kLo, kHi, ~, openHi] = certaincoast.betweensamples.timeScaleInterval(t, x, kGrid);
    [lo, hi] = certaincoast.betweensamples.crossingBand(kLo, kHi, kGrid, curve);
    note = '';
    if openHi
        note = ' (shortest time scales not excluded)';
    end
    fprintf('%-18s correlation time %6.1f d, 95%% interval %6.1f-%6.1f d%s; warning band %.2f-%.2f\n', ...
            designs{d, 1}, 1 / kHat, 1 / kHi, 1 / kLo, note, lo, hi);
end
