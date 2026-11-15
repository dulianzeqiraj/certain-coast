function arms = aquiferArms(M, pump, S)
%AQUIFERARMS Five competing structures of the synthetic aquifer model: the arms of TRIAL-GO.
%   ARMS = certaincoast.calibration.aquiferArms(M, PUMP) returns a 1-by-5 struct array for the aquifer M
%   (certaincoast.model.planAquifer) and the pumping PUMP (fields Q0 and Q, see aquiferOutputs). Each structure maps
%   its parameters theta (log-multipliers with the Gaussian prior N(0, Gam)) to a field of conductivity
%   multipliers u = B theta(1:nK) and to the recharge M.recharge exp(theta(end)):
%       A1  zones, recharge fixed          one multiplier per pilot zone, prior M.zonePrior; recharge not estimated,
%                                          as in common PEST practice
%       A2  uniform conductivity            one multiplier (standard deviation 1) and the recharge
%       A3  three facies                    one multiplier per facies class (standard deviation 0.5) and the recharge
%       A4  zones + recharge                A1 with the recharge estimated
%       A5  zones + coastal strip           A4 and one multiplier of the cells within M.stripKm of the coast
%                                          (standard deviation 1)
%   The log recharge has the standard deviation S.sigR (default 0.3). The fields of an arm are id, name, B, nK,
%   hasR, th0, Gam and fwd, the handle theta -> certaincoast.calibration.aquiferOutputs; ARMS can be passed to
%   certaincoast.calibration.trialGo directly.
%   ARMS = certaincoast.calibration.aquiferArms(M, PUMP, S) uses the settings S (certaincoast.calibration.defaults).
%
%   See also certaincoast.calibration.trialGo, certaincoast.calibration.aquiferTruth.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 3 || isempty(S)
    S = certaincoast.calibration.defaults();
else
    S = certaincoast.calibration.defaults(S);
end
Nc = M.Nc;
Bz = sparse((1:Nc).', M.zone, 1, Nc, M.Nz);
Bf = sparse((1:Nc).', M.facies, 1, Nc, 3);
Bs = sparse(find(M.strip), 1, 1, Nc, 1);
def = {'A1', 'zones, recharge fixed',    Bz,                  M.zonePrior,                false;
       'A2', 'uniform conductivity',     sparse(ones(Nc, 1)), 1,                          true;
       'A3', 'three facies',             Bf,                  0.25 * eye(3),              true;
       'A4', 'zones + recharge',         Bz,                  M.zonePrior,                true;
       'A5', 'zones + coastal strip',    [Bz Bs],             blkdiag(M.zonePrior, 1),    true};
for m = 1:size(def, 1)
    B = def{m, 3};
    nK = size(B, 2);
    hasR = def{m, 5};
    if hasR
        Gam = blkdiag(def{m, 4}, S.sigR ^ 2);
        factor = @(th) exp(th(end));
    else
        Gam = def{m, 4};
        factor = @(th) 1;
    end
    a = struct('id', def{m, 1}, 'name', def{m, 2}, 'B', B, 'nK', nK, 'hasR', hasR, ...
               'th0', zeros(size(Gam, 1), 1), 'Gam', Gam);
    a.fwd = @(th) certaincoast.calibration.aquiferOutputs(M, B * th(1:nK), M.recharge * factor(th), pump);
    arms(m) = a; %#ok<AGROW>
end
end
