%EXAMPLE_MOTION Vertical land motion and relative sea level from InSAR seen from two orbits (synthetic data).
%   A synthetic coastal point is seen by Sentinel-1 from an ascending and a descending orbit every 6 days from 2016 to
%   2026. The ground sinks steadily, swells and shrinks with the seasons, and steps down in an earthquake on
%   26 November 2019, with a slower post-seismic settling. Each line-of-sight series is fitted (rate, step, seasonal
%   cycle, post-seismic term, AR(1) errors); the two looks are combined into vertical and east-west motion; the
%   vertical rate is turned into relative sea level. No real site or data is used.
%   To run the same steps on public data, download an InSAR product (for example from the European Ground Motion
%   Service) and read it with certaincoast.motion.readSeriesCsv. Run from the repository root.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

addpath(fileparts(fileparts(mfilename('fullpath'))));
rng(7);
quake = certaincoast.motion.decimalYear(2019, 11, 26);
t = certaincoast.motion.decimalYear(2016, 1, 1) + (0:6:3900).' / 365.25;
n = numel(t);

% synthetic truth (mm): vertical and east motion
vUp = -4.0; vEast = 1.0;                                       % mm/yr
stepUp = -15; stepEast = 5;                                    % coseismic steps
post = -3; tau = 0.1;                                          % post-seismic settling, mm per unit log, years
up = vUp * (t - 2016) + stepUp * (t >= quake) + 3 * sin(2 * pi * (t - 0.2)) ...
     + post * log(1 + max(t - quake, 0) / tau);
east = vEast * (t - 2016) + stepEast * (t >= quake);

% line-of-sight unit vectors [east north up] towards the satellite (incidence about 39 degrees)
uA = [-0.61 -0.14 0.78]; uA = uA / norm(uA);
uD = [0.61 -0.14 0.78];  uD = uD / norm(uD);
ar1 = @(phi, s) filter(1, [1 -phi], s * sqrt(1 - phi ^ 2) * randn(n, 1));   % correlated noise
losA = uA(1) * east + uA(3) * up + ar1(0.5, 2.5);
losD = uD(1) * east + uD(3) * up + ar1(0.5, 2.5);

opts = struct('steps', quake, 'postseismic', struct('epoch', quake, 'tau', tau));
FA = certaincoast.motion.fitMotionSeries(t, losA, opts);
FD = certaincoast.motion.fitMotionSeries(t, losD, opts);
fprintf('ascending  LOS: rate %5.2f +- %.2f mm/yr, step %6.2f +- %.2f mm (AR(1) inflation %.1f)\n', ...
        FA.rate, FA.rateSigma, FA.step, FA.stepSigma, FA.inflation);
fprintf('descending LOS: rate %5.2f +- %.2f mm/yr, step %6.2f +- %.2f mm (AR(1) inflation %.1f)\n', ...
        FD.rate, FD.rateSigma, FD.step, FD.stepSigma, FD.inflation);

[vu, ve, su, se] = certaincoast.motion.losToVertical(FA.rate, FD.rate, uA, uD, FA.rateSigma, FD.rateSigma);
[du, de] = certaincoast.motion.losToVertical(FA.step, FD.step, uA, uD);
fprintf('vertical rate %5.2f +- %.2f mm/yr (true %.1f); east rate %5.2f +- %.2f mm/yr (true %.1f)\n', ...
        vu, su, vUp, ve, se, vEast);
fprintf('coseismic step: vertical %6.2f mm (true %d), east %5.2f mm (true %d)\n', du, stepUp, de, stepEast);

% relative sea level: a sea rising 3.3 +- 0.4 mm/yr (an assumed value for the example) over this land
[r, s, rise, riseSigma] = certaincoast.motion.relativeSeaLevel(3.3, 0.4, vu, su, 2050 - 2026);
fprintf('relative sea level %4.1f +- %.1f mm/yr; %3.0f +- %.0f mm more by 2050 at constant rates\n', ...
        r, s, rise, riseSigma);
