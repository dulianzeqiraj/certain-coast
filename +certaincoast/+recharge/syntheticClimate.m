function [P, PET] = syntheticClimate(nDays, startDay)
%SYNTHETICCLIMATE Synthetic daily rain and potential evapotranspiration for a Mediterranean coastal plain.
%   [P, PET] = certaincoast.recharge.syntheticClimate(NDAYS) returns NDAYS days of rain P (mm/d; wet winters, dry
%   summers, gamma-distributed amounts) and potential evapotranspiration PET (mm/d; 0.7 in January, 5.0 in July),
%   starting on 1 September.
%   [...] = certaincoast.recharge.syntheticClimate(NDAYS, STARTDAY) starts on day of year STARTDAY (1 = 1 January).
%
%   Synthetic data only: no real station is used. The gamma draws use only core functions.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 2
    startDay = 245;
end
doy = mod((0:nDays - 1).' + startDay - 1, 365);
s = 0.5 * (1 + cos(2 * pi * (doy - 350) / 365));          % wetness of the season, peak mid-December
pWet = 0.06 + 0.40 * s;
meanAmount = 5 + 9 * s;
wet = rand(nDays, 1) < pWet;
amount = gammaShapeBelowOne(0.75, nDays) .* meanAmount / 0.75;
P = wet .* amount;
PET = 0.7 + 4.3 * 0.5 * (1 - cos(2 * pi * (doy - 15) / 365));
end

function g = gammaShapeBelowOne(k, n)
% Gamma(k, 1) draws for 0 < k < 1: G = Gamma(k + 1) * U^(1/k), Gamma(k + 1) by Marsaglia and Tsang (2000).
d = k + 1 - 1 / 3;
c = 1 / sqrt(9 * d);
g = zeros(n, 1);
for i = 1:n
    while true
        x = randn;
        v = (1 + c * x) ^ 3;
        if v > 0 && log(rand) < 0.5 * x ^ 2 + d - d * v + d * log(v)
            break
        end
    end
    g(i) = d * v * rand ^ (1 / k);
end
end
