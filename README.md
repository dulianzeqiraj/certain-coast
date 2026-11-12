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
| `certaincoast.model` | steady variable-density flow and salt transport in a coastal section, with a sparse coloured Jacobian (Henry benchmark); compaction, liquefaction and coastal water level follow in WP2 | available (flow and salt) |
| `certaincoast.network` | monitoring and recharge design for the decision: goal-oriented information gain from two adjoint solves, greedy and exhaustive design, facies prior and pooling loss, information about a yes/no question (is a supply well salinised?), recharge sites that lower that risk | available |
| `certaincoast.surrogate` | fast machine-learned copy of the model with a certified error bound | planned (WP2) |
| `certaincoast.decisions` | safe operating window, pumping rule and safe pumping tests | planned (WP5) |

## Requirements

MATLAB R2016b or later, or GNU Octave 7 or later. **No toolbox is needed**: only core functions are used
(`betaincinv`, `fminbnd`, `fzero`, `lsqnonneg`, `erfcinv`).

## Install and test

```matlab
% from the folder of the repository
addpath(pwd)
run_tests                          % 36 tests, about 20 s
run examples/example_alarm.m
run examples/example_between_samples.m
run examples/example_boundary.m
run examples/example_recharge.m
run examples/example_henry.m
run examples/example_network.m
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

## Origin of the model and network packages

`certaincoast.model` and `certaincoast.network` come from the author's code for goal-oriented experimental design of
seawater intrusion and its application to the Fushë-Kuqe coastal aquifer (Zeqiraj, D., 2026, *Journal of Hydrology*
680, 136559, https://doi.org/10.1016/j.jhydrol.2026.136559). They are rewritten here without any field data:
the residual is vectorised (identical to the original to the last digit), the Jacobian is sparse and coloured, and
every test runs on the synthetic Henry section.

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
Environmental Science Center of Qatar University as partner and the Ministry of Infrastructure and Energy of Albania as
end-user.
