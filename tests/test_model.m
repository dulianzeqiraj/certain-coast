function tests = test_model()
%TEST_MODEL Tests of certaincoast.model (run them with run_tests).
tests = {@henry_toe_matches_the_benchmark_values, @coloured_jacobian_matches_column_differences, ...
         @salt_and_water_balance_close, @adjoint_matches_finite_differences};
end

function henry_toe_matches_the_benchmark_values()
% toe of the 0.5 isochlor at the base: 1.4226 m (16x8) and 1.3971 m (32x16), as in Zeqiraj (2026)
P = certaincoast.model.henryParameters();
expected = [1.4226, 1.3971];
sizes = [16 8; 32 16];
for k = 1:2
    g = certaincoast.model.henryGrid(sizes(k, 1), sizes(k, 2));
    [U, normF] = certaincoast.model.solveDensity(zeros(g.N, 1), g, P, 1e-12);
    toe = certaincoast.model.toePosition(U, g);
    assert(normF < 1e-10 && abs(toe - expected(k)) < 5e-5, 'toe %.4f m, residual %.1e', toe, normF);
end
end

function coloured_jacobian_matches_column_differences()
rng(12);
g = certaincoast.model.henryGrid(7, 5);
P = certaincoast.model.henryParameters();
U = [0.05 * randn(g.N, 1); rand(g.N, 1)];
theta = 0.5 * randn(g.N, 1);
[J, F0, Jt] = certaincoast.model.densityJacobian(U, theta, g, P);
Jd = zeros(2 * g.N);
for k = 1:2 * g.N
    d = 1e-7 * max(1, abs(U(k)));
    Up = U;
    Up(k) = Up(k) + d;
    Jd(:, k) = (certaincoast.model.densityResidual(Up, theta, g, P) - F0) / d;
end
assert(norm(full(J) - Jd, 'fro') <= 1e-12 * norm(Jd, 'fro'));
k = 9;
d = 1e-7 * max(1, abs(theta(k)));
tp = theta;
tp(k) = tp(k) + d;
col = (certaincoast.model.densityResidual(U, tp, g, P) - F0) / d;
assert(norm(full(Jt(:, k)) - col) <= 1e-12 * max(norm(col), 1));
end

function salt_and_water_balance_close()
% at steady state the fresh water entering inland leaves at the sea
g = certaincoast.model.henryGrid(16, 8);
P = certaincoast.model.henryParameters();
U = certaincoast.model.solveDensity(zeros(g.N, 1), g, P, 1e-12);
h = reshape(U(1:g.N), g.Nx, g.Ny).';
K = P.Kref * ones(g.Ny, 1);
qSea = -K .* (P.eps * (g.Ly - g.yc(:)) - h(:, g.Nx)) / (0.5 * g.dx);
assert(abs(sum(qSea) * g.dy - P.Q * g.Ly) < 1e-10 * P.Q);
end

function adjoint_matches_finite_differences()
% sensitivities of the salinity at three cells and of the salt mass from two adjoint solves, against central
% differences of the full nonlinear solution in a random direction
rng(13);
g = certaincoast.model.henryGrid(8, 4);
P = certaincoast.model.henryParameters();
theta = 0.3 * randn(g.N, 1);
U = certaincoast.model.solveDensity(theta, g, P, 1e-13);
[J, ~, Jt] = certaincoast.model.densityJacobian(U, theta, g, P);
obs = [5, 14, 27];
dOdu = sparse(1:3, g.N + obs, 1, 3, 2 * g.N);
dqdu = [zeros(1, g.N), g.dx * g.dy * ones(1, g.N)];
[G, H] = certaincoast.network.adjointSensitivities(J, Jt, dOdu, dqdu);
v = randn(g.N, 1);
v = v / norm(v);
delta = 1e-5;
Up = certaincoast.model.solveDensity(theta + delta * v, g, P, 1e-13);
Um = certaincoast.model.solveDensity(theta - delta * v, g, P, 1e-13);
fdObs = (Up(g.N + obs) - Um(g.N + obs)) / (2 * delta);
fdGoal = (dqdu * (Up - Um)) / (2 * delta);
assert(norm(G * v - fdObs) <= 1e-5 * norm(fdObs), 'observations: %.2e', norm(G * v - fdObs) / norm(fdObs));
assert(abs(H * v - fdGoal) <= 1e-5 * abs(fdGoal), 'goal: %.2e', abs(H * v - fdGoal) / abs(fdGoal));
end
