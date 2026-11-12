function [lo, hi, delta] = certifyBand(lo, hi, calibLo, calibHi, calibTruth, alpha)
%CERTIFYBAND Widen a robust band so that it covers the truth with a stated probability, even outside the model family.
%   [LO, HI, DELTA] = certaincoast.betweensamples.certifyBand(LO, HI, CALIBLO, CALIBHI, CALIBTRUTH, ALPHA) widens the
%   band [LO, HI] of a probability by the split-conformal margin DELTA of the distances by which the calibration
%   truths CALIBTRUTH fell outside their bands [CALIBLO, CALIBHI] (zero when inside). For a new site exchangeable with
%   the calibration sites, the widened band covers the truth with probability at least 1 - ALPHA (Proposition S3 of
%   the preparatory experiment PATH-GO). The widened band is clipped to [0, 1].
%
%   See also certaincoast.alarms.conformalMargin, certaincoast.betweensamples.crossingBand.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

d = max(max(calibLo(:) - calibTruth(:), calibTruth(:) - calibHi(:)), 0);
delta = certaincoast.alarms.conformalMargin(d, alpha);
lo = max(lo - delta, 0);
hi = min(hi + delta, 1);
end
