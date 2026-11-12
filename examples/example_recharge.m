%EXAMPLE_RECHARGE Recharge from the heads of one well, and why one anchor is needed.
%   Three years of synthetic Mediterranean rain run through the soil-water balance to the head of a well; the head
%   logger has 1 cm noise. The water-table fluctuation estimate gives daily and monthly recharge when the specific
%   yield is known; scaling five parameters together leaves the heads unchanged but changes the recharge.
%   Synthetic data only. Run from the repository root.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

addpath(fileparts(fileparts(mfilename('fullpath'))));
rng(11);
par = struct('m', 1.12, 'k_et', 0.9, 'Ic', 30, 'theta_fc', 150, 'tau_s', 8, 'phi0', 0.25, 'Sy', 0.2, ...
             'tau_a', 150, 'h_b', 0.8, 'h_0', 1.4);
N = 1096;
[P, PET] = certaincoast.recharge.syntheticClimate(N);
Pg = P / par.m;                                           % the gauge catches less than the true rain
o = certaincoast.recharge.soilWaterBalance(par, Pg, PET);
fprintf('rain %.0f mm/yr, recharge %.0f mm/yr (%.0f%% of the rain)\n', sum(P) / 3, sum(o.R) / 3, ...
        100 * sum(o.R) / sum(P));

h = o.h + 0.01 * randn(N, 1);
R = certaincoast.recharge.rechargeFromHeads(h, par.Sy, par.tau_a, par.h_b, par.h_0);
monthly = @(x) accumarray(ceil((1:numel(x)).' / 30), x(:));      % sums over 30-day blocks
fprintf('recharge from heads: daily error %.1f mm (predicted %.1f); monthly error %.1f mm (mean monthly %.1f)\n', ...
        sqrt(mean((R(2:end) - o.R(2:end)) .^ 2)), 1000 * par.Sy * 0.01 * sqrt(2), ...
        sqrt(mean((monthly(R(2:end)) - monthly(o.R(2:end))) .^ 2)), mean(monthly(o.R(2:end))));

q = par;
names = {'m', 'k_et', 'Ic', 'theta_fc', 'Sy'};
for i = 1:numel(names)
    q.(names{i}) = 1.5 * q.(names{i});
end
o2 = certaincoast.recharge.soilWaterBalance(q, Pg, PET);
fprintf('scaled parameters: largest head change %.1e m, recharge %.0f mm/yr instead of %.0f\n', ...
        max(abs(o2.h - o.h)), sum(o2.R) / 3, sum(o.R) / 3);
fprintf('so a pumping test (Sy), a soil-moisture probe or a chloride balance is needed to fix the amount.\n');
