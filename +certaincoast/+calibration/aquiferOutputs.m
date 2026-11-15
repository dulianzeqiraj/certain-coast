function [f, H] = aquiferOutputs(M, u, recharge, pump)
%AQUIFEROUTPUTS Outputs of one run of the synthetic aquifer: heads, metered pumping stops and the goal.
%   F = certaincoast.calibration.aquiferOutputs(M, U, RECHARGE, PUMP) runs certaincoast.model.planHeads once for
%   the aquifer M (certaincoast.model.planAquifer), the log-multiplier field U of the conductivity and the recharge
%   RECHARGE [m/d], and returns the column of 3J + 4 outputs, J = numel(M.obsCell):
%       F(1:J)          heads at the observation wells under the current pumping PUMP.Q0 = [QA0 QB0] [m3/d]
%                       (the calibration data of PEST practice)
%       F(J+1:2J)       head rise at the observation wells when the pumping at A is stopped (metered change QA0)
%       F(2J+1:3J)      head rise when the pumping at B is stopped (metered change QB0)
%       F(3J+1:3J+2)    future heads at A and B under the planned pumping PUMP.Q = [QA QB] (the goal)
%       F(3J+3:3J+4)    heads at A and B under recharge alone: the change of the future heads per unit change of
%                       the recharge factor
%   All pumping states are right-hand sides of the same run. A stop needs no new well, only a meter and the
%   wells already read: by reciprocity its response at well j is the sensitivity of the head at the supply well to
%   pumping at j, and by flux-scale invariance it fixes the scale of the conductivity that heads under recharge
%   cannot fix (certaincoast.model.planHeads).
%   [F, H] = ... also returns the heads [hR, gA, gB] of certaincoast.model.planHeads.
%
%   See also certaincoast.calibration.aquiferArms, certaincoast.calibration.aquiferTruth.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

H = certaincoast.model.planHeads(M, u, recharge);
o = M.obsCell(:);
s = M.supplyCell(:);
hR = H(:, 1);
gA = H(:, 2);
gB = H(:, 3);
f = [hR(o) + pump.Q0(1) * gA(o) + pump.Q0(2) * gB(o);
     -pump.Q0(1) * gA(o);
     -pump.Q0(2) * gB(o);
     hR(s) + pump.Q(1) * gA(s) + pump.Q(2) * gB(s);
     hR(s)];
end
