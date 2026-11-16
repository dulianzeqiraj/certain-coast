# Changelog

## 0.2.0 (2026-11-16)
* `certaincoast.calibration` (new): calibration when the model may be wrong (TRIAL-GO). Competing model structures
  are the arms of a bandit; a play is one Levenberg-Marquardt iteration with a prior (`lmStart`, `lmStep`,
  `lmPredict`), checked against a validation item the structure was not fitted to; the robust Gittins index of
  Cohen and Treetanthiploet (2022, Appendix D; `robustGittinsTable`, `robustGittinsIndex`) decides the next play
  (`trialGo`, with UCB, Thompson and greedy indices for comparison). Synthetic benchmark: five structures of a coastal
  aquifer, truths richer than all of them, metered pumping stops, future heads and intrusion risk at two supply
  wells (`aquiferArms`, `aquiferTruth`, `aquiferOutputs`, `aquiferGoal`).
* `certaincoast.model`: plan-view confined coastal aquifer (`planAquifer`, synthetic) and its steady heads from one
  factorisation (`planHeads`), with reciprocity and flux-scale invariance.
* `figures/`: scripts of the publication figures of the examples (`run_all_figures.m`), of the robust certified
  sensor network (`robust_certified_network.m`) and of the diagrams of TRIAL-GO and of the whole solution
  (`trial_go_scheme.m`, `solution_flowchart.m`).
* `certaincoast.motion` (new): land motion and relative sea level from InSAR and GNSS. Fit of a displacement series
  (rate, steps such as an earthquake, annual and semi-annual cycle, logarithmic post-seismic term) with AR(1)
  errors, so that rates are not over-confident (`fitMotionSeries`); vertical and east-west motion from two
  orbits (`losToVertical`); relative sea-level rate and its extrapolation (`relativeSeaLevel`); reader of public
  InSAR products with one column per date, such as EGMS (`readSeriesCsv`); `decimalYear`.
* 18 new tests (54 in all), `examples/example_trial_go.m` and `examples/example_motion.m`.

## 0.1.0 (2026-11-12)
* First public release, in MATLAB; the same code runs in GNU Octave. No toolbox is needed.
* `certaincoast.alarms`: split-conformal alarm margin for a stated miss rate; exact (Clopper-Pearson) bound on the
  miss rate observed in operation.
* `certaincoast.betweensamples`: probability of crossing a limit between samples (exact bridge maximum); time scale
  and its likelihood-ratio interval from irregular records; robust warning band; certified widening; Fisher
  information per pair of samples.
* `certaincoast.boundary`: exact passive inland boundary; reduced passive boundary by non-negative least squares;
  positive exponential sum for a long exterior; boundary learned from two piezometers; full and truncated simulators.
* `certaincoast.recharge`: soil-water balance from rain to heads; recharge from heads; synthetic climate.
* `certaincoast.model`: steady variable-density flow and salt transport in a vertical section (Henry benchmark),
  vectorised residual, sparse Jacobian by coloured finite differences (also with respect to the log-conductivity),
  damped Newton solver, toe position.
* `certaincoast.network`: goal-oriented information operator and gain, adjoint sensitivities, cost-aware greedy
  design, facies prior with a positivity check, pooling loss, information about a yes/no goal (closed form and exact
  ensemble estimate), risk-oriented recharge siting. From the author's Fushë-Kuqe code, without any field data.
* `certaincoast.util`: bounded nonlinear least squares; bivariate normal distribution function.
* 36 tests (`run_tests`), six examples, continuous integration in MATLAB and in GNU Octave.
