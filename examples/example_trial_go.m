%EXAMPLE_TRIAL_GO Calibration when the model may be wrong: TRIAL-GO, a robust bandit over five model structures.
%   Synthetic coastal aquifer of certaincoast.model.planAquifer: no real well, asset or location is used. Every
%   truth is richer than every structure: a sub-zonal conductivity field that no structure can represent and, in
%   half of the truths, a low-conductivity coastal strip that only structure A5 has. The data are the heads at 25
%   observation wells under the current pumping and the metered head changes of two pumping stops at the supply
%   wells (no drilling). The goal is the future head at the two supply wells under the planned pumping, and the
%   risk that either is intruded (head below 1 m); the decision is "safe" when that risk is below 0.1.
%   Methods compared on the same truths:
%     PEST practice      structure A1 (zones, recharge fixed) calibrated on the heads
%     PEST on all data   structure A4 (zones and recharge) calibrated on the heads and the stops
%     TRIAL-GO           robust Gittins index over A1-A5 (k = 0.5, beta = 0.95), then the chosen one on all data
%     BIC                all five structures on the heads, the one of smallest BIC, then on all data
%   About a minute in GNU Octave. Run from the repository root.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

addpath(fileparts(fileparts(mfilename('fullpath'))));
M = certaincoast.model.planAquifer();
S = certaincoast.calibration.defaults();
J = numel(M.obsCell);

% planned pumping: brings the prior-mean head at each supply well to the critical head; current pumping: half of it
H = certaincoast.model.planHeads(M, zeros(M.Nc, 1), M.recharge);
c = M.supplyCell;
Q = [H(c, 2) H(c, 3)] \ (M.hc - H(c, 1));
pump = struct('Q', Q.', 'Q0', 0.5 * Q.');
arms = certaincoast.calibration.aquiferArms(M, pump, S);
S.table = certaincoast.calibration.robustGittinsTable(S.k, S.beta, S.horizon);
rowsAll = 1:3 * J;
sigAll = [S.sigH * ones(J, 1); S.sigV * ones(2 * J, 1)];

% budget of TRIAL-GO: 2.5 x the median cost of calibrating A4 on the heads, over five pilot truths
pilot = zeros(1, 5);
for i = 1:5
    tr = certaincoast.calibration.aquiferTruth(M, pump, 1000 + i, false, S);
    st = certaincoast.calibration.lmStart(arms(4), 1:J, tr.dH, S.sigH, S);
    while ~st.conv
        st = certaincoast.calibration.lmStep(arms(4), st, S);
    end
    pilot(i) = st.runs;
end
S.budget = round(2.5 * median(pilot));

nT = 30;
names = {'PEST practice (A1, heads)', 'PEST on all data (A4)', 'TRIAL-GO, then all data', 'BIC, then all data'};
cover = zeros(nT, 2, 4);
crps = zeros(nT, 4);
wrong = zeros(nT, 4);
runs = zeros(nT, 4);
chosen = zeros(nT, 1);
strip = false(nT, 1);
missRate = zeros(nT, 5);
for t = 1:nT
    strip(t) = mod(t, 2) == 0;
    tr = certaincoast.calibration.aquiferTruth(M, pump, t, strip(t), S);
    cal = struct('rows', 1:J, 'd', tr.dH, 'sigma', S.sigH);
    val = struct('rows', J + 1:3 * J, 'd', tr.dV, 'sigma', S.sigV);
    dAll = [tr.dH; tr.dV];
    % every structure on the heads: PEST practice (A1), the BIC and the miss rate of the stops at convergence
    heads = cell(1, 5);
    bic = zeros(1, 5);
    cost = 0;
    for m = 1:5
        st = certaincoast.calibration.lmStart(arms(m), cal.rows, cal.d, cal.sigma, S);
        while ~st.conv
            st = certaincoast.calibration.lmStep(arms(m), st, S);
        end
        heads{m} = st;
        cost = cost + st.runs;
        bic(m) = sum(((tr.dH - st.f(1:J)) / S.sigH) .^ 2) + numel(arms(m).th0) * log(J);
        [mu, C] = certaincoast.calibration.lmPredict(st, val.rows);
        missRate(t, m) = mean(abs(val.d - mu) > S.z * sqrt(max(diag(C), 0) + S.sigV ^ 2));
    end
    goal = cell(1, 4);
    goal{1} = certaincoast.calibration.aquiferGoal(heads{1}, M, S);
    runs(t, 1) = heads{1}.runs;
    st = certaincoast.calibration.lmStart(arms(4), rowsAll, dAll, sigAll, S);
    while ~st.conv
        st = certaincoast.calibration.lmStep(arms(4), st, S);
    end
    goal{2} = certaincoast.calibration.aquiferGoal(st, M, S);
    runs(t, 2) = st.runs;
    rng(100000 + t);                                   % the order in which each structure sees the stops
    R = certaincoast.calibration.trialGo(arms, cal, val, S);
    goal{3} = certaincoast.calibration.aquiferGoal(R.final, M, S);
    runs(t, 3) = R.runsTotal;
    chosen(t) = R.sel;
    [~, mb] = min(bic);
    st = certaincoast.calibration.lmStart(arms(mb), rowsAll, dAll, sigAll, S, heads{mb});
    while ~st.conv
        st = certaincoast.calibration.lmStep(arms(mb), st, S);
    end
    goal{4} = certaincoast.calibration.aquiferGoal(st, M, S);
    runs(t, 4) = cost + st.runs;
    for q = 1:4
        g = goal{q};
        z = (tr.Hfut - g.mu) ./ g.s;
        cover(t, :, q) = (abs(z) <= S.z).';
        crps(t, q) = mean(g.s .* (z .* erf(z / sqrt(2)) + 2 * exp(-z .^ 2 / 2) / sqrt(2 * pi) - 1 / sqrt(pi)));
        wrong(t, q) = (g.risk < 0.1) ~= (tr.riskTrue < 0.1);
    end
end

fprintf('%d truths (%d with the coastal strip); budget of TRIAL-GO %d runs (2.5 x %d)\n', nT, sum(strip), ...
        S.budget, median(pilot));
fprintf('%-27s %9s %9s %10s %16s %6s\n', 'method', 'cover A', 'cover B', 'CRPS (m)', 'wrong decisions', 'runs');
for q = 1:4
    fprintf('%-27s %9.2f %9.2f %10.3f %16.2f %6.0f\n', names{q}, mean(cover(:, 1, q)), mean(cover(:, 2, q)), ...
            mean(crps(:, q)), mean(wrong(:, q)), mean(runs(:, q)));
end
fprintf('miss rate of the stops after calibration on the heads (nominal 0.10), A1-A5:\n');
fprintf('  without the strip %s\n  with the strip    %s\n', mat2str(mean(missRate(~strip, :)), 2), ...
        mat2str(mean(missRate(strip, :)), 2));
fprintf('structures chosen by TRIAL-GO, A1-A5: without the strip %s, with the strip %s\n', ...
        mat2str(accumarray(chosen(~strip), 1, [5 1]).'), mat2str(accumarray(chosen(strip), 1, [5 1]).'));
