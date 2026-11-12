function tests = test_alarms()
%TEST_ALARMS Tests of certaincoast.alarms (run them with run_tests).
tests = {@margin_gives_the_stated_coverage, @too_few_errors_give_an_infinite_margin, @margin_is_an_order_statistic, ...
         @alarm_rule, @clopper_pearson_bound_with_no_miss, @clopper_pearson_bound_with_one_miss, @invalid_input};
end

function margin_gives_the_stated_coverage()
rng(0);
m = certaincoast.alarms.conformalMargin(randn(2000, 1), 0.10);
coverage = mean(randn(20000, 1) <= m);
assert(coverage >= 0.89, 'coverage %.3f below the stated 0.90', coverage);
end

function too_few_errors_give_an_infinite_margin()
assert(isinf(certaincoast.alarms.conformalMargin([0.1 0.2 0.3], 0.05)));
end

function margin_is_an_order_statistic()
% n = 19, alpha = 0.1: k = ceil(20 * 0.9) = 18, the 18th smallest error
r = (19:-1:1) / 10;
assert(abs(certaincoast.alarms.conformalMargin(r, 0.10) - 1.8) < 1e-12);
end

function alarm_rule()
assert(isequal(certaincoast.alarms.issueAlarm([1.0 2.0], 2.5, 0.6), [false true]));
end

function clopper_pearson_bound_with_no_miss()
u = certaincoast.alarms.missRateUpperBound(0, 10, 0.95);
assert(abs(u - (1 - 0.05 ^ (1 / 10))) < 1e-9);
end

function clopper_pearson_bound_with_one_miss()
% the bound solves P(Binomial(20, u) <= 1) = 0.05
u = certaincoast.alarms.missRateUpperBound(1, 20);
cdf = (1 - u) ^ 20 + 20 * u * (1 - u) ^ 19;
assert(abs(cdf - 0.05) < 1e-9);
end

function invalid_input()
assert(throws(@() certaincoast.alarms.conformalMargin([], 0.1)));
assert(throws(@() certaincoast.alarms.conformalMargin([1 2 3], 1.5)));
assert(throws(@() certaincoast.alarms.missRateUpperBound(3, 2)));
end

function tf = throws(f)
tf = false;
try
    f();
catch
    tf = true;
end
end
