function tests = test_motion()
%TEST_MOTION Tests of certaincoast.motion (run them with run_tests).
tests = {@decimal_years, @fit_recovers_rate_step_and_season, @correlated_errors_widen_the_rate, ...
         @los_decomposition_is_exact, @relative_sea_level_adds_subsidence, @reads_a_product_with_date_columns, ...
         @invalid_input};
end

function decimal_years()
assert(abs(certaincoast.motion.decimalYear(2020, 1, 1) - (2020 + 0.5 / 366)) < 1e-12);
assert(abs(certaincoast.motion.decimalYear(2021, 12, 31) - (2021 + 364.5 / 365)) < 1e-12);
assert(certaincoast.motion.decimalYear(20191126) == certaincoast.motion.decimalYear(2019, 11, 26));
end

function fit_recovers_rate_step_and_season()
rng(1);
t = 2016 + (0:6:3650).' / 365.25;
d = 2 - 4 * (t - 2016) - 15 * (t >= 2019.9) + 3 * cos(2 * pi * t) + 1.5 * randn(size(t));
F = certaincoast.motion.fitMotionSeries(t, d, struct('steps', 2019.9, 'noise', 'white'));
assert(abs(F.rate + 4) < 3 * F.rateSigma && F.rateSigma < 0.2, 'rate %.3f +- %.3f', F.rate, F.rateSigma);
assert(abs(F.step + 15) < 3 * F.stepSigma, 'step %.2f +- %.2f', F.step, F.stepSigma);
assert(abs(F.annualAmplitude - 3) < 0.6, 'annual amplitude %.2f', F.annualAmplitude);
assert(abs(F.sigma0 - 1.5) < 0.2 && F.inflation == 1);
end

function correlated_errors_widen_the_rate()
rng(2);
t = 2016 + (0:6:3650).' / 365.25;
e = filter(1, [1 -0.8], 2 * randn(size(t)));
d = -4 * (t - 2016) + e;
W = certaincoast.motion.fitMotionSeries(t, d, struct('noise', 'white', 'seasonal', false));
R = certaincoast.motion.fitMotionSeries(t, d, struct('seasonal', false));
assert(R.phi > 0.6 && R.inflation > 4, 'phi %.2f', R.phi);
assert(R.rateSigma > 2 * W.rateSigma, 'AR(1) did not widen the rate');
end

function los_decomposition_is_exact()
uA = [-0.61 -0.14 0.78]; uA = uA / norm(uA);
uD = [0.61 -0.14 0.78]; uD = uD / norm(uD);
vE = [1 -2 0.5]; vU = [-4 3 0];
losA = uA(1) * vE + uA(3) * vU;
losD = uD(1) * vE + uD(3) * vU;
[u, e, su, se] = certaincoast.motion.losToVertical(losA, losD, uA, uD, 1, 1);
assert(max(abs(u - vU)) < 1e-12 && max(abs(e - vE)) < 1e-12);
assert(all(su > 0) && all(se > 0) && isequal(size(su), size(vU)));
end

function relative_sea_level_adds_subsidence()
[r, s, rise, riseSigma] = certaincoast.motion.relativeSeaLevel(3.3, 0.4, -4.0, 0.3, 10);
assert(abs(r - 7.3) < 1e-12 && abs(s - 0.5) < 1e-12);
assert(abs(rise - 73) < 1e-9 && abs(riseSigma - 5) < 1e-9);
end

function reads_a_product_with_date_columns()
f = [tempname() '.csv'];
fid = fopen(f, 'w');
fprintf(fid, 'pid,easting,northing,mean_velocity,20200101,20200113,20200125\n');
fprintf(fid, 'A1,100.5,200.25,-3.1,0.0,-0.4,-1.1\n');
fprintf(fid, 'B2,101.5,201.25,1.2,0.0,0.2,0.5\n');
fclose(fid);
S = certaincoast.motion.readSeriesCsv(f);
delete(f);
assert(isequal(size(S.D), [2 3]) && abs(S.D(1, 3) + 1.1) < 1e-12);
assert(abs(S.t(1) - certaincoast.motion.decimalYear(2020, 1, 1)) < 1e-12 && isequal(S.dates(2), 20200113));
assert(iscell(S.attr.pid) && strcmp(S.attr.pid{2}, 'B2'));
assert(isnumeric(S.attr.easting) && abs(S.attr.mean_velocity(1) + 3.1) < 1e-12);
end

function invalid_input()
u = [0.61 -0.14 0.78]; u = u / norm(u);
assert(throws(@() certaincoast.motion.losToVertical(1, 2, u, u)));
assert(throws(@() certaincoast.motion.fitMotionSeries([1 2 3], [1 2 3])));
assert(throws(@() certaincoast.motion.fitMotionSeries(2016 + (0:99) / 50, randn(1, 100), struct('steps', 2030))));
assert(throws(@() certaincoast.motion.decimalYear(2020, 13, 1)));
end

function tf = throws(f)
tf = false;
try
    f();
catch
    tf = true;
end
end
