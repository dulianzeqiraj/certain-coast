function tests = test_recharge()
%TEST_RECHARGE Tests of certaincoast.recharge (run them with run_tests).
tests = {@water_balance_closes, @recharge_from_exact_heads_is_exact, @daily_error_matches_the_prediction, ...
         @heads_do_not_fix_the_amount_of_recharge, @least_squares_respects_the_bounds};
end

function par = truth()
par = struct('m', 1.12, 'k_et', 0.9, 'Ic', 30, 'theta_fc', 150, 'tau_s', 8, 'phi0', 0.25, 'Sy', 0.2, ...
             'tau_a', 150, 'h_b', 0.8, 'h_0', 1.4);
end

function water_balance_closes()
rng(7);
par = truth();
[P, PET] = certaincoast.recharge.syntheticClimate(1096);
o = certaincoast.recharge.soilWaterBalance(par, P / par.m, PET);
storeChange = o.theta(end) - par.phi0 * par.theta_fc;
assert(abs(sum(P) - sum(o.runoff) - sum(o.ET) - sum(o.R) - storeChange) < 1e-8);
end

function recharge_from_exact_heads_is_exact()
rng(8);
par = truth();
[P, PET] = certaincoast.recharge.syntheticClimate(400);
o = certaincoast.recharge.soilWaterBalance(par, P / par.m, PET);
R = certaincoast.recharge.rechargeFromHeads(o.h, par.Sy, par.tau_a, par.h_b, par.h_0);
assert(max(abs(R - o.R)) < 1e-9);
R0 = certaincoast.recharge.rechargeFromHeads(o.h, par.Sy, par.tau_a, par.h_b);
assert(isnan(R0(1)) && max(abs(R0(2:end) - o.R(2:end))) < 1e-9);
end

function daily_error_matches_the_prediction()
rng(9);
par = truth();
[P, PET] = certaincoast.recharge.syntheticClimate(1096);
o = certaincoast.recharge.soilWaterBalance(par, P / par.m, PET);
sigma = 0.01;
R = certaincoast.recharge.rechargeFromHeads(o.h + sigma * randn(1096, 1), par.Sy, par.tau_a, par.h_b, par.h_0);
e = R(2:end) - o.R(2:end);
predicted = 1000 * par.Sy * sigma * sqrt(2);
assert(abs(sqrt(mean(e .^ 2)) / predicted - 1) < 0.1);
end

function heads_do_not_fix_the_amount_of_recharge()
rng(10);
par = truth();
[P, PET] = certaincoast.recharge.syntheticClimate(730);
Pg = P / par.m;
o = certaincoast.recharge.soilWaterBalance(par, Pg, PET);
q = par;
names = {'m', 'k_et', 'Ic', 'theta_fc', 'Sy'};
for i = 1:numel(names)
    q.(names{i}) = 1.7 * q.(names{i});
end
o2 = certaincoast.recharge.soilWaterBalance(q, Pg, PET);
assert(max(abs(o2.h - o.h)) < 1e-12);                     % same heads ...
assert(abs(sum(o2.R) / sum(o.R) - 1.7) < 1e-9);           % ... with 70% more recharge
end

function least_squares_respects_the_bounds()
% minimise (x1 - 2)^2 + (x2 + 1)^2 with x >= 0: the answer is (2, 0)
fun = @(X) [X(1, :) - 2; X(2, :) + 1];
x = certaincoast.util.boundedLeastSquares(fun, [1; 1], [0; 0]);
assert(norm(x - [2; 0]) < 1e-6);
end
