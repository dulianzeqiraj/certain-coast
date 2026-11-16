% certaincoast.motion
%   Land motion and relative sea level from InSAR and GNSS: the vertical motion of the coast drives seawater
%   intrusion, the water table under quays and foundations, and the reach of storm surges at the same time.
%
%   fitMotionSeries  - rate, steps (earthquakes), seasonal cycle and post-seismic term of a displacement series,
%                      with AR(1) errors so that the rate is not over-confident
%   losToVertical    - vertical and east-west motion from line-of-sight motion seen from two orbits
%   relativeSeaLevel - relative sea-level rate (sea rise minus vertical land motion) and its extrapolation
%   readSeriesCsv    - read public InSAR products with one point per row and one column per date (EGMS, LiCSBAS)
%   decimalYear      - dates as decimal years
%
% Start with examples/example_motion.m (synthetic data). The repository holds no field data: public products are
% downloaded by the user, and results for individual assets are not published here (DATA_POLICY.md).
%
% Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.
