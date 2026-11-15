function S = defaults(S)
%DEFAULTS Settings of certaincoast.calibration: TRIAL-GO, Levenberg-Marquardt and the synthetic aquifer.
%   S = certaincoast.calibration.defaults() returns the default settings.
%   S = certaincoast.calibration.defaults(S) adds the default of every field missing from S.
%
%   Levenberg-Marquardt (lmStart, lmStep)
%       fdStep     1e-3     forward-difference step in the parameters (log-multipliers)
%       lam0       1e-2     initial damping; divided by 4 after an accepted step, multiplied by 8 after a rejected one
%       tolRel     1e-4     converged when an accepted step lowers the objective by less than this fraction
%       maxIt      15       converged after this many iterations (also when the damping exceeds 1e6)
%   Bandit (trialGo)
%       budget     Inf      model runs the bandit may spend
%       policy     'trial'  index: 'trial' (robust Gittins), 'ucb', 'thompson' or 'greedy'
%       k          0.5      credible level of the set of miss rates of the robust index (k = 0: no ambiguity)
%       beta       0.95     discount of the robust index
%       horizon    60       plays looked ahead by the robust index
%       table      []       robust Gittins table (robustGittinsTable); computed from k, beta, horizon when empty
%       z          1.6449   half-width of the central interval in standard deviations (90%) for the check of an item
%       minItems   5        validation items an arm needs before it can be selected
%       order      []       order in which every arm sees the validation items (one row per arm); random when empty
%       condition  true     calibrate the selected arm on all data at the end
%   Synthetic aquifer (aquiferTruth, aquiferGoal)
%       sigH       0.05     noise of the heads [m]
%       sigV       0.01     noise of the metered head changes of a pumping stop [m]
%       sigR       0.3      standard deviation of the log recharge (truth and prior of the arms)
%       sigXi      0.08     standard deviation of the future change of the recharge factor
%       zetaSd     0.35     standard deviation of the sub-zonal log-conductivity field of a truth
%       zetaL      2        correlation length of that field [km] (exponential covariance)
%       stripLog   -1       log-multiplier of the low-conductivity coastal strip of a truth
%
%   See also certaincoast.calibration.trialGo, certaincoast.calibration.lmStart.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

d = struct();
d.fdStep = 1e-3;
d.lam0 = 1e-2;
d.tolRel = 1e-4;
d.maxIt = 15;
d.budget = Inf;
d.policy = 'trial';
d.k = 0.5;
d.beta = 0.95;
d.horizon = 60;
d.table = [];
d.z = 1.6448536269514722;
d.minItems = 5;
d.order = [];
d.condition = true;
d.sigH = 0.05;
d.sigV = 0.01;
d.sigR = 0.3;
d.sigXi = 0.08;
d.zetaSd = 0.35;
d.zetaL = 2;
d.stripLog = -1;
if nargin < 1 || isempty(S)
    S = d;
    return
end
if ~isstruct(S) || ~isscalar(S)
    error('certaincoast:calibration:badSettings', 'S must be a scalar structure.');
end
names = fieldnames(d);
for i = 1:numel(names)
    if ~isfield(S, names{i})
        S.(names{i}) = d.(names{i});
    end
end
end
