function R = trialGo(arms, cal, val, S)
%TRIALGO Calibration when the model may be wrong: a robust bandit over competing model structures (TRIAL-GO).
%   R = certaincoast.calibration.trialGo(ARMS, CAL, VAL, S) calibrates competing structures of a model (the arms of
%   a bandit) on the data CAL and selects one, while it checks each structure against validation items VAL it has
%   not been fitted to and spends a fixed budget of model runs S.budget where it is best spent.
%   ARMS is a struct array, one element per structure, with the fields
%       fwd     handle: the column of all outputs of one model run for the parameters theta, f = fwd(theta)
%       th0     prior mean of theta;   Gam   prior covariance of theta (also the Tikhonov regularisation)
%       name    (optional) a label
%   CAL has the fields rows (indices of the outputs that are data), d (values) and sigma (noise standard deviation,
%   a vector or a scalar); VAL the same for the validation items, one item per value, for example the metered head
%   changes of a pumping stop. S holds the settings (certaincoast.calibration.defaults).
%
%   A play of arm m is one Levenberg-Marquardt iteration of m on CAL (certaincoast.calibration.lmStep: a finite-
%   difference Jacobian of n_m runs and one trial run; the first play also evaluates the prior mean; once m has
%   converged a play is free) followed by the check of the next item of VAL that m has not seen. The cost of the
%   play is a miss,
%       h = 1  if  |y - mu| > z sqrt(v + sigma^2),   otherwise h = 0,
%   with mu and v the mean and variance of the linearised predictive of the item at the current iterate of m
%   (certaincoast.calibration.lmPredict) and z = S.z (1.645: the central 90% interval). After n plays with s misses
%   the state of m is its posterior-mean miss rate p = (s + 1/2)/(n + 1). Every arm is played once; then the arm of
%   smallest index is played, as long as some arm has items left and its next play fits in the budget. The index of
%   S.policy = 'trial' is the robust Gittins index of Cohen and Treetanthiploet (2022), read from S.table at
%   (p, n + 1) (certaincoast.calibration.robustGittinsIndex): the miss rate with an exploration bonus at small n and
%   a premium for the ambiguity of a short record. For comparison, 'ucb' uses the lower confidence bound
%   p - sqrt(2 log N / n) (N plays in all), 'thompson' a draw from Beta(s + 1/2, n - s + 1/2) and 'greedy' p.
%   When no play fits any more, the arm of smallest p among the arms with at least S.minItems items is selected
%   (ties: smaller index, then smaller arm number); with S.condition true (the default) it is then calibrated on
%   all data, CAL and VAL, from its current iterate, and these runs are added in R.runsTotal.
%
%   R has the fields
%       sel, name        selected arm and its label;   policy   the index used
%       n, s, p          items checked, misses and posterior-mean miss rate of every arm
%       runs, runsTotal  model runs of the bandit, and with the final calibration on all data
%       states           calibration state of every arm at the end of the bandit (lmStart, lmStep)
%       final            state of the selected arm calibrated on all data ([] when S.condition is false)
%       log              one row per play: [arm, runs so far, miss];  logIndex  the index of every arm before it
%       order            order in which every arm saw the items (S.order; by default each arm gets its own random
%                        permutation from the current random stream, so set it with rng to repeat a run)
%
%   Example (two structures of a straight line; the second one is right)
%       x = linspace(0, 1, 40).';  y = 1 + 2 * x + 0.05 * randn(40, 1);
%       arms(1) = struct('name', 'constant', 'fwd', @(t) t(1) + 0 * x, 'th0', 0, 'Gam', 10);
%       arms(2) = struct('name', 'line', 'fwd', @(t) t(1) + t(2) * x, 'th0', [0; 0], 'Gam', 10 * eye(2));
%       cal = struct('rows', 1:2:40, 'd', y(1:2:40), 'sigma', 0.05);
%       val = struct('rows', 2:2:40, 'd', y(2:2:40), 'sigma', 0.05);
%       R = certaincoast.calibration.trialGo(arms, cal, val, struct('budget', 60));
%
%   Reference: S. N. Cohen and T. Treetanthiploet (2022), Gittins' theorem under uncertainty, Electronic Journal of
%   Probability 27, paper 17, https://doi.org/10.1214/22-EJP742
%
%   See also certaincoast.calibration.robustGittinsTable, certaincoast.calibration.aquiferArms.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 4 || isempty(S)
    S = certaincoast.calibration.defaults();
else
    S = certaincoast.calibration.defaults(S);
end
nA = numel(arms);
if nA < 1 || ~isstruct(arms) || ~all(isfield(arms, {'fwd', 'th0', 'Gam'}))
    error('certaincoast:calibration:badArms', 'ARMS must be a struct array with the fields fwd, th0 and Gam.');
end
cal = dataSet(cal, 'CAL');
val = dataSet(val, 'VAL');
nItems = numel(val.d);
policy = lower(S.policy);
if ~any(strcmp(policy, {'trial', 'ucb', 'thompson', 'greedy'}))
    error('certaincoast:calibration:badPolicy', 'unknown policy %s.', S.policy);
end
order = S.order;
if isempty(order)
    order = zeros(nA, nItems);
    for m = 1:nA
        order(m, :) = randperm(nItems);
    end
elseif ~isequal(size(order), [nA nItems]) || any(any(sort(order, 2) ~= repmat(1:nItems, nA, 1)))
    error('certaincoast:calibration:badOrder', 'S.order must hold one permutation of 1..%d per arm.', nItems);
end
Tab = S.table;
if strcmp(policy, 'trial') && isempty(Tab)
    Tab = certaincoast.calibration.robustGittinsTable(S.k, S.beta, S.horizon);
end

nPar = zeros(1, nA);
for m = 1:nA
    nPar(m) = numel(arms(m).th0);
end
st = cell(1, nA);
started = false(1, nA);
n = zeros(1, nA);
s = zeros(1, nA);
used = 0;
logPlay = zeros(0, 3);
logIndex = zeros(0, nA);
first = true;
while true
    ok = false(1, nA);
    for m = 1:nA
        ok(m) = n(m) < nItems && playCost(st{m}, started(m), nPar(m)) <= S.budget - used;
    end
    if first                                       % every arm once, in order
        first = false;
        g = nan(1, nA);
        for m = 1:nA
            if n(m) < nItems && playCost(st{m}, started(m), nPar(m)) <= S.budget - used
                [st{m}, h, c] = playArm(arms(m), st{m}, started(m), cal, val, order(m, n(m) + 1), S);
                started(m) = true;
                used = used + c;
                n(m) = n(m) + 1;
                s(m) = s(m) + h;
                logPlay(end + 1, :) = [m used h]; %#ok<AGROW>
                logIndex(end + 1, :) = g; %#ok<AGROW>
            end
        end
        continue
    end
    if ~any(ok)
        break
    end
    g = inf(1, nA);
    for m = find(ok)
        g(m) = armIndex(policy, Tab, s(m), n(m), sum(n));
    end
    [~, m] = min(g);
    [st{m}, h, c] = playArm(arms(m), st{m}, started(m), cal, val, order(m, n(m) + 1), S);
    started(m) = true;
    used = used + c;
    n(m) = n(m) + 1;
    s(m) = s(m) + h;
    logPlay(end + 1, :) = [m used h]; %#ok<AGROW>
    logIndex(end + 1, :) = g; %#ok<AGROW>
end

% selection: smallest posterior-mean miss rate among the arms with enough items
p = (s + 0.5) ./ (n + 1);
gf = inf(1, nA);
if ~strcmp(policy, 'thompson')
    for m = find(n > 0)
        gf(m) = armIndex(policy, Tab, s(m), n(m), sum(n));
    end
end
elig = find(n >= S.minItems);
if isempty(elig)
    [~, sel] = max(n);
else
    [~, o] = sortrows([p(elig).' gf(elig).' elig(:)]);
    sel = elig(o(1));
end

% the selected arm on all data, from its current iterate
final = [];
runsFinal = 0;
if S.condition
    rowsAll = [cal.rows; val.rows];
    dAll = [cal.d; val.d];
    sigAll = [cal.sigma; val.sigma];
    if started(sel)
        final = certaincoast.calibration.lmStart(arms(sel), rowsAll, dAll, sigAll, S, st{sel});
    else
        final = certaincoast.calibration.lmStart(arms(sel), rowsAll, dAll, sigAll, S);
    end
    while ~final.conv
        final = certaincoast.calibration.lmStep(arms(sel), final, S);
    end
    runsFinal = final.runs;
end
if isfield(arms, 'name') && ~isempty(arms(sel).name)
    name = arms(sel).name;
else
    name = sprintf('arm %d', sel);
end
R = struct('sel', sel, 'name', name, 'policy', policy, 'n', n, 's', s, 'p', p, 'runs', used, ...
           'runsTotal', used + runsFinal);
R.states = st;
R.final = final;
R.log = logPlay;
R.logIndex = logIndex;
R.order = order;
end

function D = dataSet(D, label)
% check a data set {rows, d, sigma} and return it as columns
if ~isstruct(D) || ~all(isfield(D, {'rows', 'd', 'sigma'}))
    error('certaincoast:calibration:badData', '%s must have the fields rows, d and sigma.', label);
end
D.rows = D.rows(:);
D.d = D.d(:);
D.sigma = D.sigma(:);
if isscalar(D.sigma)
    D.sigma = D.sigma + zeros(size(D.d));
end
if isempty(D.d) || numel(D.rows) ~= numel(D.d) || numel(D.sigma) ~= numel(D.d) || any(~(D.sigma > 0))
    error('certaincoast:calibration:badData', '%s needs as many rows and positive sigmas as values.', label);
end
end

function c = playCost(st, started, nPar)
% model runs of the next play of an arm
if ~started
    c = nPar + 2;                                  % prior mean, Jacobian, trial run
elseif st.conv
    c = 0;
else
    c = isempty(st.J) * nPar + 1;
end
end

function [st, h, c] = playArm(arm, st, started, cal, val, q, S)
% one Levenberg-Marquardt iteration of the arm, then the check of validation item q
c = 0;
if ~started
    st = certaincoast.calibration.lmStart(arm, cal.rows, cal.d, cal.sigma, S);
    c = st.runs;
end
st = certaincoast.calibration.lmStep(arm, st, S);
c = c + st.costLast;
[mu, v] = certaincoast.calibration.lmPredict(st, val.rows(q));
h = double(abs(val.d(q) - mu) > S.z * sqrt(max(v, 0) + val.sigma(q) ^ 2));
end

function g = armIndex(policy, Tab, s, n, nAll)
% index of an arm after n plays with s misses (the smallest index is played)
p = (s + 0.5) / (n + 1);
switch policy
    case 'trial'
        g = certaincoast.calibration.robustGittinsIndex(Tab, p, n + 1);
    case 'ucb'
        g = p - sqrt(2 * log(max(nAll, 1)) / max(n, 1));
    case 'thompson'
        g = betaDraw(s + 0.5, n - s + 0.5, rand);
    otherwise
        g = p;
end
end

function x = betaDraw(a, b, u)
% Beta(a, b) draw by inversion of the uniform u (bisection on betainc; no toolbox needed)
lo = 0;
hi = 1;
for it = 1:50
    x = 0.5 * (lo + hi);
    if betainc(x, a, b) < u
        lo = x;
    else
        hi = x;
    end
end
x = 0.5 * (lo + hi);
end
