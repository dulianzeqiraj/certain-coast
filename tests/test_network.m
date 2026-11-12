function tests = test_network()
%TEST_NETWORK Tests of certaincoast.network and certaincoast.util.bvnCdf (run them with run_tests).
tests = {@gain_equals_mutual_information, @greedy_picks_the_sensor_that_sees_the_goal, @facies_prior_limits, ...
         @pooling_costs_nothing_when_facies_are_one, @bivariate_normal_against_closed_form, ...
         @recharge_risk_matches_finite_differences, @threshold_gain_matches_the_ensemble_estimate};
end

function gain_equals_mutual_information()
% Lemma 1: for a linear-Gaussian model the gain is 1/2 log(prior variance / posterior variance) of the goal
rng(21);
n = 6;
m = 3;
A = randn(n);
Gam = A * A' + n * eye(n);
G = randn(m, n);
H = randn(1, n);
Gn = 0.5 * eye(m);
post = inv(inv(Gam) + G' * (Gn \ G));
mi = 0.5 * log((H * Gam * H') / (H * post * H'));
assert(abs(certaincoast.network.eigOperator(G, H, Gam, Gn) - mi) < 1e-10);
end

function greedy_picks_the_sensor_that_sees_the_goal()
n = 5;
Gcand = eye(n);
H = [0 0 1 0 0];                               % the goal is parameter 3
S = certaincoast.network.greedyDesign(Gcand, H, eye(n), 0.1, 1);
assert(isequal(S, 3));
end

function facies_prior_limits()
% correlation blocks, as in Zeqiraj (2026): positive for every rhoF, independent facies at 0; invalid blocks refused
blocks = {0.6 * eye(3) + 0.4 * ones(3), 0.5 * eye(2) + 0.5 * ones(2), 0.7 * eye(4) + 0.3 * ones(4)};
for rho = [0 0.5 0.99]
    [~, flag] = chol(certaincoast.network.faciesPrior(rho, blocks));
    assert(flag == 0);
end
assert(min(eig(certaincoast.network.faciesPrior(1, blocks))) > -1e-12);
assert(isequal(certaincoast.network.faciesPrior(0, blocks), blkdiag(blocks{:})));
refused = false;
try
    certaincoast.network.faciesPrior(1, {eye(3), 0.2 * eye(4)});
catch
    refused = true;
end
assert(refused);
end

function pooling_costs_nothing_when_facies_are_one()
rng(23);
blocks = {eye(3), eye(3)};
Gcand = randn(6, 6);
H = [1 1 1 -1 -1 -1] / 3;
out = certaincoast.network.poolingLoss(Gcand, H, blocks, [], 0.2, 2, [0 0.5 1]);
assert(all(out.L >= 0) && out.L(end) < 1e-12);
end

function bivariate_normal_against_closed_form()
for r = [-0.9 -0.3 0 0.4 0.95]
    assert(abs(certaincoast.util.bvnCdf(0, 0, r) - (0.25 + asin(r) / (2 * pi))) < 1e-10);
end
Phi = @(x) 0.5 * erfc(-x / sqrt(2));
assert(abs(certaincoast.util.bvnCdf(0.7, -1.2, 0) - Phi(0.7) * Phi(-1.2)) < 1e-10);
end

function recharge_risk_matches_finite_differences()
% one and two supply wells on a small confined aquifer: closed-form gradient against finite differences
n = 30;
A = spdiags([-ones(n, 1), 2.2 * ones(n, 1), -ones(n, 1)], -1:1, n, n);
mu = [1.1; 0.9];
C = [0.09 0.03; 0.03 0.04];
hc = 1.0;
cells = [8; 21];
[risk, info] = certaincoast.network.rechargeRisk(A, cells, C, mu, hc);
P = @(m) 1 - certaincoast.util.bvnCdf((m(1) - hc) / sqrt(C(1, 1)), (m(2) - hc) / sqrt(C(2, 2)), ...
                                       C(1, 2) / sqrt(C(1, 1) * C(2, 2)));
dQ = 1e-6;
for k = [3 8 15 21 28]
    fd = -(P(mu + dQ * info.g(k, :).') - P(mu - dQ * info.g(k, :).')) / (2 * dQ);
    assert(abs(risk(k) - fd) < 1e-6 * max(abs(fd), 1e-3));
end
[risk1, info1] = certaincoast.network.rechargeRisk(A, cells(1), C(1, 1), mu(1), hc);
P1 = @(m) 0.5 * erfc((m - hc) / sqrt(2 * C(1, 1)));
fd = -(P1(mu(1) + dQ * info1.g(12)) - P1(mu(1) - dQ * info1.g(12))) / (2 * dQ);
assert(abs(risk1(12) - fd) < 1e-6 * max(abs(fd), 1e-3));
end

function threshold_gain_matches_the_ensemble_estimate()
% closed form (quadrature) against likelihood re-weighting of a Gaussian ensemble, one and two goal wells
rng(24);
m = 4;
B = randn(m + 2);
Sigma = B * B' / (m + 2) + 0.2 * eye(m + 2);
mu = [0.3; -0.2];
goal = [m + 1, m + 2];
design = [1 3];
s = 0.3;
L = chol(Sigma, 'lower');
Zin = (L * randn(m + 2, 30000)).';
Zout = (L * randn(m + 2, 3000)).';
Zin(:, goal) = bsxfun(@plus, Zin(:, goal), mu.');
Yout = Zout(:, 1:m) + s * randn(3000, m);
gains = [Zin(:, goal(1)) < 0, (Zin(:, goal(1)) < 0) | (Zin(:, goal(2)) < 0)];
ex = certaincoast.network.exactGoalEig(Zin(:, 1:m), gains, Yout, design, s);
one = certaincoast.network.thresholdGoalEig(Sigma, mu(1), design, s ^ 2, goal(1), 0);
two = certaincoast.network.thresholdGoalEig(Sigma, mu, design, s ^ 2, goal, 0);
assert(abs(one - ex(1)) < 0.02 && abs(two - ex(2)) < 0.02, 'closed %.3f %.3f, ensemble %.3f %.3f', one, two, ex(1), ex(2));
end
