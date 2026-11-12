function tests = test_boundary()
%TEST_BOUNDARY Tests of certaincoast.boundary (run them with run_tests).
tests = {@exact_exterior_reproduces_the_whole_aquifer, @exterior_is_passive, @reduced_boundary_beats_plain_ones, ...
         @energy_never_increases, @half_order_sum_reaches_its_tolerance, @learned_boundary_predicts_a_new_regime};
end

function [p, H, h0] = setup()
p = certaincoast.boundary.aquifer();
A = certaincoast.boundary.conductanceMatrix(p.nFull, p, true);
f = zeros(p.nFull, 1);
f(1:p.nModel) = p.R0 * p.dx;
f(p.well) = f(p.well) - p.Q0;
hs = A \ f;
H = hs(p.boundaryCell);
h0 = hs(1:p.nModel);
end

function exact_exterior_reproduces_the_whole_aquifer()
[p, H, h0] = setup();
[~, r, Q] = certaincoast.boundary.forcing(p, 3 * 365);
full = certaincoast.boundary.simulateFull(p, r, Q);
[C, w, mu] = certaincoast.boundary.exactExterior(p);
h = certaincoast.boundary.simulateTruncated(p, r, Q, C, w, mu, H, h0);
assert(max(max(abs(h - full(:, 1:p.nModel)))) < 1e-9);
end

function exterior_is_passive()
p = certaincoast.boundary.aquifer();
[C, w, mu, Craw] = certaincoast.boundary.exactExterior(p);
assert(C >= 0 && all(w >= 0) && all(mu > 0) && abs(Craw) < 1e-9);  % closed far end: no steady flux
end

function reduced_boundary_beats_plain_ones()
[p, H, h0] = setup();
[~, r, Q] = certaincoast.boundary.forcing(p, 5 * 365);
full = certaincoast.boundary.simulateFull(p, r, Q);
[C, w, mu] = certaincoast.boundary.exactExterior(p);
[Cf, wf, muf] = certaincoast.boundary.fitPositive(@(s) certaincoast.boundary.impedance(s, C, w, mu), 12, ...
                                                 [2 * pi / (40 * 365.25), 2 * pi / 5]);
bc = p.boundaryCell;
err = @(h) sqrt(mean((h(:, bc) - full(:, bc)) .^ 2));
e12 = err(certaincoast.boundary.simulateTruncated(p, r, Q, Cf, wf, muf, H, h0));
eNoFlow = err(certaincoast.boundary.simulateTruncated(p, r, Q, 0, 0, 1, H, h0));
eFixed = err(certaincoast.boundary.simulateTruncated(p, r, Q, p.c, 0, 1, H, h0));
assert(all(wf >= 0) && Cf >= 0);
assert(e12 < 0.01 * min(eNoFlow, eFixed), 'errors %.4f, %.4f, %.4f', e12, eNoFlow, eFixed);
end

function energy_never_increases()
p = certaincoast.boundary.aquifer();
[C, w, mu] = certaincoast.boundary.exactExterior(p);
n = 2 * 365;
h0 = zeros(p.nModel, 1);
h0(9:13) = 1;
pp = p;
pp.R0 = 0;
[~, ~, E] = certaincoast.boundary.simulateTruncated(pp, zeros(n, 1), zeros(n, 1), C, w, mu, 0, h0);
assert(all(diff(E) <= 1e-12));
end

function half_order_sum_reaches_its_tolerance()
wmin = 2 * pi / (30 * 365.25);
[w, mu, M] = certaincoast.boundary.sqrtImpedanceSum(wmin, wmin * 1e4, 1e-2);
om = logspace(log10(wmin), log10(wmin * 1e4), 2000).';
e = max(abs(certaincoast.boundary.impedance(1i * om, 0, w, mu) - sqrt(1i * om)) ./ abs(sqrt(1i * om)));
assert(M == 14 && e < 0.02, 'M = %d, error %.4f', M, e);
end

function learned_boundary_predicts_a_new_regime()
rng(6);
[p, H, h0] = setup();
nd = 8 * 365;
[~, r, Q] = certaincoast.boundary.forcing(p, nd, 'broadband');
full = certaincoast.boundary.simulateFull(p, r, Q);
o1 = full(:, p.bandCells(1)) + 0.002 * randn(nd, 1);
o2 = full(:, p.bandCells(2)) + 0.002 * randn(nd, 1);
[C, w, mu] = certaincoast.boundary.identifyFromPiezometers(p, o1, o2, r, H);
assert(C >= 0 && all(w >= 0));
% a drought (recharge x 0.3 in years 3-8) and twice the pumping from year 2, never seen in the calibration
[~, rp, Qp] = certaincoast.boundary.forcing(p, 10 * 365, 'seasonal', [3 * 365, 8 * 365, 0.3], [2 * 365, Inf, 1.2]);
ref = certaincoast.boundary.simulateFull(p, rp, Qp);
bc = p.boundaryCell;
err = @(h) sqrt(mean((h(:, bc) - ref(:, bc)) .^ 2));
eLearned = err(certaincoast.boundary.simulateTruncated(p, rp, Qp, C, w, mu, H, h0));
eNoFlow = err(certaincoast.boundary.simulateTruncated(p, rp, Qp, 0, 0, 1, H, h0));
assert(eLearned < 0.05 * eNoFlow, 'learned %.4f m, no-flow %.4f m', eLearned, eNoFlow);
end
