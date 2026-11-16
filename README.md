# CERTAIN-COAST

**Certified early warning for coastal aquifers and infrastructure.**
Open-source MATLAB code of the CERTAIN-COAST project (the coupled coastal model and the alarm algorithm;
a fast machine-learned copy of the model follows in a later version), written by Dulian Zeqiraj (Faculty of Geology and Mining, Polytechnic University of Tirana).
The same code runs in the free GNU Octave, so a water utility or a ministry can use it without a MATLAB licence.

[![tests](https://github.com/dulianzeqiraj/certain-coast/actions/workflows/tests.yml/badge.svg)](https://github.com/dulianzeqiraj/certain-coast/actions)
<!-- After the first release on Zenodo, add the DOI badge here. -->

## What it does

On a low, subsiding coast, seawater intrusion into drinking-water wells, liquefaction in an earthquake and extreme
sea levels are linked through the relative sea level and the water table. CERTAIN-COAST watches them together and
issues warnings that **state in advance how often they may miss a dangerous event**, a rate that is then checked
against every logged outcome.

| Package | Content | Status |
|---|---|---|
| `certaincoast.alarms` | alarm margin calibrated for a stated miss rate; exact bound on the miss rate observed in operation | available |
| `certaincoast.betweensamples` | probability that chloride crosses the limit *between* laboratory samples; robust band over all models the record allows; certified widening; value of a logger or of paired samples | available |
| `certaincoast.boundary` | inland boundary of the model replaced by a passive Robin term and storage compartments (exact, reduced, or learned from two piezometers) | available |
| `certaincoast.recharge` | soil-water balance from rain to the head of a well; recharge from heads; the anchor needed to fix its amount | available |
| `certaincoast.model` | steady variable-density flow and salt transport in a coastal section, with a sparse coloured Jacobian (Henry benchmark); a synthetic plan-view confined coastal aquifer for calibration and design tests; compaction, liquefaction and coastal water level follow in WP2 | available (flow and salt) |
| `certaincoast.network` | monitoring and recharge design for the decision: goal-oriented information gain from two adjoint solves, greedy and exhaustive design, facies prior and pooling loss, information about a yes/no question (is a supply well salinised?), recharge sites that lower that risk | available |
| `certaincoast.calibration` | calibration when the model may be wrong (TRIAL-GO): competing model structures are the arms of a bandit; each play is one Levenberg-Marquardt iteration, checked against a metered pumping stop the structure was not fitted to; the robust Gittins index of Cohen and Treetanthiploet (2022) spends the budget of model runs on the structures that predict | available |
| `certaincoast.motion` | land motion and relative sea level from InSAR and GNSS: rate, earthquake steps, seasons and post-seismic settling of a displacement series with honest (AR(1)) errors; vertical and east-west motion from two orbits; relative sea-level rate; reader of public InSAR products (EGMS, LiCSBAS) | available |
| `certaincoast.surrogate` | fast machine-learned copy of the model with a certified error bound | planned (WP2) |
| `certaincoast.decisions` | safe operating window, pumping rule and safe pumping tests | planned (WP5) |

## Requirements

MATLAB R2016b or later, or GNU Octave 7 or later. **No toolbox is needed**: only core functions are used
(`betaincinv`, `betainc`, `fminbnd`, `fzero`, `lsqnonneg`, `erfcinv`, `interp2`, `datenum`).

## Install and test

```matlab
% from the folder of the repository
addpath(pwd)
run_tests                          % 54 tests, about 30 s
run examples/example_alarm.m
run examples/example_between_samples.m
run examples/example_boundary.m
run examples/example_recharge.m
run examples/example_henry.m
run examples/example_network.m
run examples/example_trial_go.m      % about a minute
run examples/example_motion.m
```

In GNU Octave, from a terminal: `octave --eval "run_tests"`.

## Examples

```matlab
% an alarm with a 10% miss rate, checked in operation
m     = certaincoast.alarms.conformalMargin(pastObserved - pastForecast, 0.10);
alarm = certaincoast.alarms.issueAlarm(newForecast, 250, m);         % chloride, mg/L
u     = certaincoast.alarms.missRateUpperBound(1, 20);               % 95% upper bound after 1 miss in 20 events

% does chloride cross 250 mg/L between two samples?
kGrid = exp(linspace(log(1/400), log(4), 41));
curve = certaincoast.betweensamples.crossingProbability(kGrid, (250 - 180)/35, 120);
[~, kLo, kHi] = certaincoast.betweensamples.timeScaleInterval(sampleDays, standardisedChloride, kGrid);
[pLo, pHi] = certaincoast.betweensamples.crossingBand(kLo, kHi, kGrid, curve);

% calibration when the model may be wrong: five structures of a synthetic aquifer and a budget of model runs
M    = certaincoast.model.planAquifer();
arms = certaincoast.calibration.aquiferArms(M, pump);                % pump: current and planned pumping (m3/d)
cal  = struct('rows', 1:25,  'd', heads, 'sigma', 0.05);             % heads under the current pumping
val  = struct('rows', 26:75, 'd', stops, 'sigma', 0.01);             % metered head changes of two pumping stops
R    = certaincoast.calibration.trialGo(arms, cal, val, struct('budget', 213));
g    = certaincoast.calibration.aquiferGoal(R.final, M);              % future heads at the supply wells; risk

% land motion from InSAR seen from two orbits, and relative sea level
FA = certaincoast.motion.fitMotionSeries(t, losAsc, struct('steps', quake));   % rate, step, seasons, AR(1) errors
FD = certaincoast.motion.fitMotionSeries(t, losDesc, struct('steps', quake));
[vUp, ~, sUp] = certaincoast.motion.losToVertical(FA.rate, FD.rate, uAsc, uDesc, FA.rateSigma, FD.rateSigma);
[rsl, sRsl] = certaincoast.motion.relativeSeaLevel(seaRate, seaSigma, vUp, sUp);
```

Output of the examples (synthetic data):

* **Alarm:** margin 23.9 mg/L; 23 events, none missed; 95% upper bound on the miss rate 0.12.
* **Between samples:** with a correlation time of 8 days, chloride crosses the limit within a season with
  probability 0.79, while the samples alone suggest 0.023. Warning band: quarterly samples 0.18-1.00, monthly
  0.44-1.00, weekly 0.78-0.89, daily logger 0.79-0.82.
* **Boundary:** over 20 years, error at the boundary cell 1e-13 m with the exact passive boundary, 0.0003 m with
  12 compartments, 0.23 m with no flow and 0.25 m with a fixed head. Learned from two piezometers, the boundary
  predicts a drought with doubled pumping within 0.005 m (no flow: 0.61 m; fixed head: 1.23 m).
* **Henry benchmark:** toe of the 0.5 isochlor at 1.4226, 1.3971, 1.3918 and 1.3864 m on grids of 16 x 8, 32 x 16,
  41 x 21 and 61 x 31 cells (residual 3.5e-14 at 41 x 21), as in Zeqiraj (2026); the coloured sparse Jacobian
  needs 27 residual evaluations instead of 2N + 1, so 41 x 21 takes half a second.
* **Network:** three of twelve salinity sensors in the Henry section with an uncertain conductivity field: the
  greedy design is the exhaustive optimum for the salt mass near the toe (2.54 nats; median design 1.46). The
  question "is a supply cell above its limit?" prefers other sensors; in this small section they carry almost the
  same information (0.499 against 0.500 nats).
* **Recharge:** from heads with 1 cm noise, daily error 2.8 mm and monthly error 2.6 mm (mean monthly recharge
  32 mm); heads alone cannot fix the amount, so one anchor (pumping test, soil probe or chloride balance) is needed.
* **Calibration (TRIAL-GO):** 30 synthetic truths of a coastal aquifer, each richer than all five structures (a
  sub-zonal field that none can represent; a low-conductivity coastal strip in 15 of them that only one has). After
  calibration on the heads, the metered stops fall outside the 90% intervals of the structures A1-A5 at the rates
  0.42, 0.78, 0.66, 0.31 and 0.26 (no strip) and 0.63, 0.91, 0.71, 0.48 and 0.38 (strip): the stops reveal the wrong
  structures that the heads let pass. The 90% intervals of the future head at the two supply wells cover the truth in
  0.47 and 0.43 of the truths for PEST practice (zones, recharge fixed, heads only), 0.63 and 0.63 for PEST on all
  data and 0.67 and 0.63 for TRIAL-GO (CRPS 0.93, 0.51 and 0.42 m; wrong safe/unsafe decisions 0.23, 0.13 and 0.13),
  with a budget of 213 model runs (266 with the final calibration on all data). All five structures with BIC
  selection, then all data, cost 427 runs and do best (CRPS 0.28 m, wrong decisions 0.03). No method reaches 0.90,
  because no structure holds the sub-zonal field. Numbers from GNU Octave 8.4 (the random streams of MATLAB differ).
* **Land motion:** a synthetic coastal point seen every 6 days from two orbits (2016-2026), sinking 4.0 mm/yr with a
  15 mm step in an earthquake, seasonal swelling and post-seismic settling: the two fitted series give a vertical
  rate of -4.13 +- 0.16 mm/yr and a vertical step of -16 mm; with correlated errors the honest rate uncertainty is
  about 1.7 times the white-noise one. With a sea rising 3.3 mm/yr (assumed), the relative sea level rises
  7.4 +- 0.4 mm/yr, 178 mm more by 2050 at constant rates. Public products (for example EGMS) are read with
  `readSeriesCsv`; results for individual assets are not published here.

## Figures

The folder `figures` holds the MATLAB scripts of the figures and diagrams of the project (synthetic data only). Each
writes PNG and vector PDF files to `results/figures`, which git does not track.

| Script | Figure |
|---|---|
| `run_all_figures.m` | the publication figures of the alarm, between-samples, boundary, recharge, Henry and network examples |
| `robust_certified_network.m` | salinity sensors placed robustly, with a certificate of optimality, when sensors share errors (common error of a data-logger hub, model error correlated over distance; sealed reference channels) |
| `trial_go_scheme.m` | diagram: how TRIAL-GO calibrates when the model may be wrong |
| `solution_flowchart.m` | flowchart: the problems of the project and how it observes, models and calibrates, warns, decides and delivers |

```matlab
addpath figures
run_all_figures                    % about 40 s in GNU Octave
robust_certified_network           % about 2.5 minutes
trial_go_scheme
solution_flowchart
```

In GNU Octave without a display: `xvfb-run -a octave --no-gui --eval "addpath figures; run_all_figures"`.

## Origin of the model, network and calibration packages

`certaincoast.model` and `certaincoast.network` come from the author's code for goal-oriented experimental design of
seawater intrusion and its application to the Fushë-Kuqe coastal aquifer (Zeqiraj, D., 2026, *Journal of Hydrology*
680, 136559, https://doi.org/10.1016/j.jhydrol.2026.136559). They are rewritten here without any field data:
the residual is vectorised (identical to the original to the last digit), the Jacobian is sparse and coloured, and
every test runs on the synthetic Henry section. `certaincoast.calibration` comes from the author's calibration
study of the same Fushë-Kuqe model (TRIAL-GO); here it runs on a synthetic plan-view aquifer
(`certaincoast.model.planAquifer`), again without any field data.

## Data and security

The repository contains code, documentation and **synthetic** test data only. Results for individual assets,
sensor locations and data provided by end-users are restricted and are never published here; see
[DATA_POLICY.md](DATA_POLICY.md).

## How to cite

Use the "Cite this repository" button on GitHub or [CITATION.cff](CITATION.cff); each release has a DOI on Zenodo.

## License

Apache License 2.0 (see [LICENSE](LICENSE) and [NOTICE](NOTICE)).

## Acknowledgement

Prepared for the NATO Science for Peace and Security (SPS) Programme application CERTAIN-COAST (call 2027-1), with the
Environmental Science Center of Qatar University as partner and the Institute of Geosciences (IGJEO) of the Polytechnic
University of Tirana as end-user.
