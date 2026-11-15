function tests = test_calibration()
%TEST_CALIBRATION Tests of certaincoast.calibration and of the plan-view aquifer (run them with run_tests).
tests = {@index_tends_to_the_miss_rate, @ambiguity_raises_the_index, @patience_lowers_the_index, ...
         @index_lookup_reads_the_table, @aquifer_reciprocity, @flux_scale_invariance, ...
         @arms_reproduce_a_truth_they_can_represent, @levenberg_marquardt_fits_a_well_specified_model, ...
         @bandit_selects_the_right_structure, @bandit_respects_the_budget, @invalid_input};
end

function index_tends_to_the_miss_rate()
% gamma - p -> 0 like 1/sqrt(n), with an exploration bonus (gamma < p) at small n
T = certaincoast.calibration.robustGittinsTable(0.5, 0.95, 60, 101, [1 20 1000], 201);
d = interp1(T.p, T.gam, 0.3) - 0.3;
assert(d(1) < -0.05, 'no exploration bonus at n = 1 (%.3f)', d(1));
assert(abs(d(3)) < 0.01 && abs(d(3)) < abs(d(2)) && abs(d(2)) < abs(d(1)), 'gamma - p does not tend to 0');
end

function ambiguity_raises_the_index()
ks = [0 0.2 0.5 0.8];
d = zeros(size(ks));
for i = 1:numel(ks)
    T = certaincoast.calibration.robustGittinsTable(ks(i), 0.95, 60, 101, 5, 201);    % n = 1..5
    d(i) = interp1(T.p, T.gam(:, end), 0.3) - 0.3;
end
assert(all(diff(d) > 0), 'the index does not grow with the ambiguity k: %s', mat2str(d, 3));
end

function patience_lowers_the_index()
bs = [0.90 0.95 0.99];
d = zeros(size(bs));
for i = 1:numel(bs)
    T = certaincoast.calibration.robustGittinsTable(0.5, bs(i), 60, 101, 5, 201);
    d(i) = interp1(T.p, T.gam(:, end), 0.3) - 0.3;
end
assert(all(diff(d) < 0), 'the index does not fall with the discount beta: %s', mat2str(d, 3));
end

function index_lookup_reads_the_table()
T = certaincoast.calibration.robustGittinsTable(0.5, 0.95, 20, 51, [2 4], 101);
g = certaincoast.calibration.robustGittinsIndex(T, T.p(7), [2 3 4 9]);
assert(abs(g(1) - T.gam(7, 1)) < 1e-12 && abs(g(3) - T.gam(7, 2)) < 1e-12);
assert(abs(g(2) - 0.5 * (T.gam(7, 1) + T.gam(7, 2))) < 1e-12, 'no linear interpolation in n');
assert(abs(g(4) - T.gam(7, 2)) < 1e-12, 'n is not clipped to the table');
end

function aquifer_reciprocity()
% the head change at well j per unit extraction at A equals the head change at A per unit extraction at j
M = certaincoast.model.planAquifer();
u = 0.4 * sin(M.X) .* cos(0.7 * M.Y);
o = M.obsCell;
Q = sparse(o, 1:numel(o), 1, M.Nc, numel(o));
H = certaincoast.model.planHeads(M, u, M.recharge, Q);
a = H(o, 2);
b = H(M.supplyCell(1), 4:end).';
assert(max(abs(a - b)) <= 1e-10 * max(abs(a)), 'reciprocity fails');
assert(all(a < 0), 'extraction must lower the heads');
end

function flux_scale_invariance()
% conductivity and recharge times c: heads under recharge unchanged, responses to pumping divided by c
M = certaincoast.model.planAquifer();
u = 0.3 * cos(M.X + M.Y);
c = 1.5;
H1 = certaincoast.model.planHeads(M, u, M.recharge);
H2 = certaincoast.model.planHeads(M, u + log(c), c * M.recharge);
assert(max(abs(H2(:, 1) - H1(:, 1))) <= 1e-9 * max(abs(H1(:, 1))), 'heads under recharge changed');
assert(max(max(abs(c * H2(:, 2:3) - H1(:, 2:3)))) <= 1e-9 * max(max(abs(H1(:, 2:3)))), 'responses do not scale');
end

function arms_reproduce_a_truth_they_can_represent()
% without the sub-zonal field, A4 (and A5 with the strip) represent the truth exactly
[M, pump, S, arms] = smallCase();
S.zetaSd = 0;
tr = certaincoast.calibration.aquiferTruth(M, pump, 3, false, S);
f = arms(4).fwd([tr.thz; tr.xiR]);
assert(max(abs(f - tr.f)) < 1e-10, 'A4 does not reproduce its own truth');
tr = certaincoast.calibration.aquiferTruth(M, pump, 4, true, S);
f = arms(5).fwd([tr.thz; S.stripLog; tr.xiR]);
assert(max(abs(f - tr.f)) < 1e-10, 'A5 does not reproduce a truth with the strip');
assert(tr.riskTrue >= 0 && tr.riskTrue <= 1);
end

function levenberg_marquardt_fits_a_well_specified_model()
[M, pump, S, arms] = smallCase();
S.zetaSd = 0;
tr = certaincoast.calibration.aquiferTruth(M, pump, 5, false, S);
J = numel(M.obsCell);
st = certaincoast.calibration.lmStart(arms(4), 1:3 * J, [tr.dH; tr.dV], [S.sigH * ones(J, 1); S.sigV * ones(2 * J, 1)], S);
while ~st.conv
    st = certaincoast.calibration.lmStep(arms(4), st, S);
end
assert(st.Phi / (3 * J) < 2, 'misfit %.2f per datum for a well-specified structure', st.Phi / (3 * J));
g = certaincoast.calibration.aquiferGoal(st, M, S);
assert(all(abs(tr.Hfut0 - g.mu) < 3 * g.s), 'goal outside three standard deviations');
assert(g.risk >= 0 && g.risk <= 1);
end

function bandit_selects_the_right_structure()
[arms, cal, val, T] = lineCase();
rng(1);
R = certaincoast.calibration.trialGo(arms, cal, val, struct('budget', 200, 'table', T));
assert(R.sel == 3, 'TRIAL-GO chose %s', R.name);
assert(R.p(3) < min(R.p(1:2)), 'the parabola does not have the smallest miss rate');
assert(R.runs <= 200 && R.runsTotal >= R.runs && ~isempty(R.final));
assert(size(R.log, 1) == sum(R.n) && sum(R.log(:, 3)) == sum(R.s));
end

function bandit_respects_the_budget()
[arms, cal, val, T] = lineCase();
pols = {'trial', 'ucb', 'thompson', 'greedy'};
for i = 1:numel(pols)
    rng(2);
    R = certaincoast.calibration.trialGo(arms, cal, val, ...
                                         struct('budget', 25, 'policy', pols{i}, 'condition', false, 'table', T));
    assert(R.runs <= 25, '%s spent %d runs', pols{i}, R.runs);
    assert(isempty(R.final) && R.runsTotal == R.runs && any(R.sel == 1:3));
end
end

function invalid_input()
assert(throws(@() certaincoast.calibration.robustGittinsTable(1.2, 0.95, 10)));
assert(throws(@() certaincoast.calibration.robustGittinsTable(0.5, 1.5, 10)));
assert(throws(@() certaincoast.model.planAquifer(struct('depth', 1))));
[arms, cal, val] = lineCase();
assert(throws(@() certaincoast.calibration.trialGo(struct([]), cal, val)));
assert(throws(@() certaincoast.calibration.trialGo(arms, cal, val, struct('policy', 'random'))));
assert(throws(@() certaincoast.calibration.lmStart(arms(1), 1:3, [1 2], 0.1)));
end

% ---------------------------------------------------------------------------------------------------- helpers
function [M, pump, S, arms] = smallCase()
M = certaincoast.model.planAquifer(struct('Lx', 6, 'Ly', 8, 'dx', 0.5, 'nObs', 12, 'supply', [2.25 5.75; 2.25 2.25]));
S = certaincoast.calibration.defaults();
H = certaincoast.model.planHeads(M, zeros(M.Nc, 1), M.recharge);
c = M.supplyCell;
Q = [H(c, 2) H(c, 3)] \ (M.hc - H(c, 1));
pump = struct('Q', Q.', 'Q0', 0.5 * Q.');
arms = certaincoast.calibration.aquiferArms(M, pump, S);
end

function [arms, cal, val, T] = lineCase()
% a constant, a line and a parabola for data from a parabola, and a small table of the robust index
rng(0);
x = linspace(0, 1, 60).';
y = 1 + 2 * x + 1.5 * x .^ 2 + 0.05 * randn(60, 1);
arms(1) = struct('name', 'constant', 'fwd', @(t) t(1) + 0 * x, 'th0', 0, 'Gam', 100);
arms(2) = struct('name', 'line', 'fwd', @(t) t(1) + t(2) * x, 'th0', [0; 0], 'Gam', 100 * eye(2));
arms(3) = struct('name', 'parabola', 'fwd', @(t) t(1) + t(2) * x + t(3) * x .^ 2, 'th0', [0; 0; 0], ...
                 'Gam', 100 * eye(3));
cal = struct('rows', 1:2:60, 'd', y(1:2:60), 'sigma', 0.05);
val = struct('rows', 2:2:60, 'd', y(2:2:60), 'sigma', 0.05);
T = certaincoast.calibration.robustGittinsTable(0.5, 0.95, 20, 51, 40, 101);
end

function tf = throws(f)
tf = false;
try
    f();
catch
    tf = true;
end
end
