function [rate, sigma, rise, riseSigma] = relativeSeaLevel(seaRate, seaSigma, vlmRate, vlmSigma, years)
%RELATIVESEALEVEL Relative sea-level rate at a site: the rise of the sea minus the vertical motion of the land.
%   [RATE, SIGMA] = certaincoast.motion.relativeSeaLevel(SEARATE, SEASIGMA, VLMRATE, VLMSIGMA) returns
%       RATE = SEARATE - VLMRATE,   SIGMA = sqrt(SEASIGMA^2 + VLMSIGMA^2)
%   in mm/yr, with SEARATE the rise of the sea surface (from altimetry, a tide gauge corrected for its own land
%   motion, or a scenario) and VLMRATE the vertical land motion, positive up (from InSAR or GNSS: subsidence is
%   negative and adds to the relative rise). The two are taken as independent.
%   [RATE, SIGMA, RISE, RISESIGMA] = certaincoast.motion.relativeSeaLevel(..., YEARS) also returns the relative rise
%   [mm] after YEARS years at constant rates (a linear extrapolation, without acceleration).
%   The relative sea level drives the three hazards of CERTAIN-COAST at once: seawater intrusion, the water table
%   under quays and foundations, and the reach of storm surges.
%
%   Example (a sea rising 3.3 mm/yr over land sinking 4 mm/yr)
%       [r, s] = certaincoast.motion.relativeSeaLevel(3.3, 0.4, -4.0, 1.0);     % 7.3 +- 1.1 mm/yr
%
%   See also certaincoast.motion.losToVertical, certaincoast.motion.fitMotionSeries.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if any(seaSigma(:) < 0) || any(vlmSigma(:) < 0)
    error('certaincoast:motion:sigma', 'standard deviations must be nonnegative.');
end
rate = seaRate - vlmRate;
sigma = sqrt(seaSigma .^ 2 + vlmSigma .^ 2);
if nargin >= 5
    rise = rate .* years;
    riseSigma = sigma .* abs(years);
end
end
