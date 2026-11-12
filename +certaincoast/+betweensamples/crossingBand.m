function [pLo, pHi] = crossingBand(kLo, kHi, kGrid, curve)
%CROSSINGBAND Band of the crossing probability over all time scales consistent with the record.
%   [PLO, PHI] = certaincoast.betweensamples.crossingBand(KLO, KHI, KGRID, CURVE) returns the smallest and largest
%   crossing probability over the rates KAPPA in [KLO, KHI], where CURVE(i) is the crossing probability at KGRID(i)
%   (from crossingProbability). The curve is interpolated in log KAPPA. The band is the robust warning: it holds for
%   every model the samples cannot rule out.
%
%   See also certaincoast.betweensamples.timeScaleInterval, certaincoast.betweensamples.certifyBand.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

lk = log(kGrid(:));
c = curve(:);
a = min(max(log(kLo), lk(1)), lk(end));
b = min(max(log(kHi), lk(1)), lk(end));
inside = lk > a & lk < b;
vals = [interp1(lk, c, a); interp1(lk, c, b); c(inside)];
pLo = min(vals);
pHi = max(vals);
end
