# Changelog

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
