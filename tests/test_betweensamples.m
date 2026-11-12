function tests = test_betweensamples()
%TEST_BETWEENSAMPLES Tests of certaincoast.betweensamples (run them with run_tests).
tests = {@bridge_maximum_counts_crossings_between_grid_points, @slow_indicator_gives_the_marginal_tail, ...
         @crossing_probability_grows_with_the_rate, @pair_information_short_lag_limit, ...
         @pair_information_matches_the_score_variance, @logger_narrows_the_band, @certified_band_widens_and_clips};
end

function bridge_maximum_counts_crossings_between_grid_points()
rng(1);
n = 2000;
X = certaincoast.betweensamples.ouPaths(1 / 8, n, 480, 0.25);
onGrid = mean(max(X, [], 2) >= 2);
withBridge = mean(certaincoast.betweensamples.bridgeMax(X, 2 / 8, 0.25) >= 2);
assert(withBridge > onGrid);
end

function slow_indicator_gives_the_marginal_tail()
rng(2);
p = certaincoast.betweensamples.crossingProbability(1e-5, 2.0, 120, 0.25, 20000);
tail = 0.5 * erfc(2 / sqrt(2));                           % 0.0228
assert(abs(p - tail) < 0.006, 'p = %.4f, marginal tail %.4f', p, tail);
end

function crossing_probability_grows_with_the_rate()
rng(3);
k = exp(linspace(log(1 / 400), log(4), 15));
P = certaincoast.betweensamples.crossingProbability(k, 2.0, 120, 0.25, 3000);
assert(all(diff(P) >= -0.01));                            % Slepian ordering, up to Monte Carlo noise
assert(P(end) > 0.95 && P(1) < 0.1);
end

function pair_information_short_lag_limit()
k = 1 / 8;
I = certaincoast.betweensamples.pairInformation(k, 1e-4) * k ^ 2;
assert(abs(I - 0.5) < 1e-3);
end

function pair_information_matches_the_score_variance()
% the variance of the score of one transition equals the Fisher information
rng(4);
k = 1 / 8;
tau = 7;
n = 200000;
x0 = randn(n, 1);
r = exp(-k * tau);
x1 = r * x0 + sqrt(1 - r ^ 2) * randn(n, 1);
h = 1e-6;
ll = @(kk) logTransition(kk, tau, x0, x1);
score = (ll(k + h) - ll(k - h)) / (2 * h);
I = certaincoast.betweensamples.pairInformation(k, tau);
assert(abs(var(score) / I - 1) < 0.03);
end

function ll = logTransition(k, tau, x0, x1)
r = exp(-k * tau);
v = 1 - r ^ 2;
ll = -0.5 * log(2 * pi * v) - (x1 - r * x0) .^ 2 / (2 * v);
end

function logger_narrows_the_band()
rng(5);
kGrid = exp(linspace(log(1 / 400), log(4), 31));
curve = certaincoast.betweensamples.crossingProbability(kGrid, 2.0, 120, 0.25, 3000);
truth = interp1(log(kGrid), curve, log(1 / 8));
width = zeros(1, 2);
steps = [91 1];
for d = 1:2
    t = certaincoast.betweensamples.sampleTimes(8 * 365, steps(d));
    x = certaincoast.betweensamples.ouAtTimes(t, 1 / 8);
    [~, kLo, kHi] = certaincoast.betweensamples.timeScaleInterval(t, x, kGrid);
    [lo, hi] = certaincoast.betweensamples.crossingBand(kLo, kHi, kGrid, curve);
    width(d) = hi - lo;
end
assert(width(2) < 0.15 && width(1) > 2 * width(2), 'widths %.2f (quarterly) and %.2f (logger)', width(1), width(2));
assert(truth > 0.6);
end

function certified_band_widens_and_clips()
[lo, hi, delta] = certaincoast.betweensamples.certifyBand(0.5, 0.95, [0.4 0.5 0.6], [0.6 0.7 0.8], ...
                                                         [0.5 0.75 0.7], 0.25);
assert(abs(delta - 0.05) < 1e-12);                        % distances 0, 0.05, 0; k = ceil(4 * 0.75) = 3
assert(abs(lo - 0.45) < 1e-12 && hi == 1);
end
