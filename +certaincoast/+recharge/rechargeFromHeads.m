function R = rechargeFromHeads(h, Sy, tauA, hB, h0)
%RECHARGEFROMHEADS Daily recharge from the heads of a well (water-table fluctuation in stochastic form).
%   R = certaincoast.recharge.rechargeFromHeads(H, SY, TAUA, HB) returns the daily recharge (mm/d) implied by the
%   daily heads H (m) for a known specific yield SY, drainage time TAUA (d) and drainage base HB (m):
%
%       R(i) = 1000 SY (H(i) - H(i-1) + (H(i-1) - HB) / TAUA).
%
%   The first value is NaN unless the head before the record, H0, is given:
%   R = certaincoast.recharge.rechargeFromHeads(H, SY, TAUA, HB, H0).
%
%   With head noise SIGMA (m), the daily error has standard deviation about 1000 SY SIGMA sqrt(2), but the errors
%   telescope: monthly and seasonal sums are much more accurate than the daily values (Proposition 1 of the
%   preparatory experiment INFILTRATION-GO).
%
%   See also certaincoast.recharge.soilWaterBalance.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

h = h(:);
if nargin >= 5 && ~isempty(h0)
    hPrev = [h0; h(1:end - 1)];
else
    hPrev = [NaN; h(1:end - 1)];
end
R = 1000 * Sy * (h - hPrev + (hPrev - hB) / tauA);
end
