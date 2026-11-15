function g = aquiferGoal(st, M, S)
%AQUIFERGOAL Goal of a calibrated structure: future heads at the supply wells and the risk of intrusion.
%   G = certaincoast.calibration.aquiferGoal(ST, M) uses the calibration ST of a structure of aquiferArms (lmStart,
%   lmStep, or R.final of trialGo) to predict the future heads at the supply wells A and B under the planned
%   pumping: the linearised predictive of lmPredict plus the variance S.sigXi^2 hnat hnat' of the future change of
%   the recharge factor (the heads are affine in it). With a = (mu - M.hc) ./ s, the risk that A or B is intruded
%   (future head below the critical head M.hc) is
%       risk = 1 - Phi2(a_A, a_B; r),
%   Phi2 the bivariate normal distribution function (certaincoast.util.bvnCdf) and r the correlation of the two
%   heads. G has the fields mu, s, r, C, lo and hi (central interval of S.z standard deviations) and risk.
%   G = certaincoast.calibration.aquiferGoal(ST, M, S) uses the settings S (certaincoast.calibration.defaults).
%
%   See also certaincoast.calibration.aquiferTruth, certaincoast.calibration.trialGo.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 3 || isempty(S)
    S = certaincoast.calibration.defaults();
else
    S = certaincoast.calibration.defaults(S);
end
J = numel(M.obsCell);
[mu, C] = certaincoast.calibration.lmPredict(st, [3 * J + 1; 3 * J + 2]);
hn = st.f(3 * J + 3:3 * J + 4);
C = C + S.sigXi ^ 2 * (hn * hn.');
C = 0.5 * (C + C.');
s = sqrt(diag(C));
r = min(max(C(1, 2) / (s(1) * s(2)), -0.999999), 0.999999);
a = (mu - M.hc) ./ s;
g = struct('mu', mu, 's', s, 'r', r, 'lo', mu - S.z * s, 'hi', mu + S.z * s, ...
           'risk', 1 - certaincoast.util.bvnCdf(a(1), a(2), r));
g.C = C;
end
