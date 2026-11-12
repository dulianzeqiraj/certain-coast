% CERTAIN-COAST: certified early warning for coastal aquifers and infrastructure.
% Version 0.1.0 (2026-11-12). MATLAB R2016b or later, or GNU Octave 7 or later; no toolbox is needed.
%
% Packages
%   certaincoast.alarms          - alarms whose miss rate is stated in advance and checked in operation
%   certaincoast.betweensamples  - warnings for limit crossings between laboratory samples (robust, certified band)
%   certaincoast.boundary        - passive inland boundary of the model, learned from two piezometers
%   certaincoast.recharge        - recharge from rain and from the heads of a well
%   certaincoast.model           - variable-density flow and salt transport in a coastal section (Henry benchmark)
%   certaincoast.network         - monitoring and recharge design for the decision: goal-oriented information gain,
%                                  adjoint sensitivities, facies prior, risk of a salinised supply well
%   certaincoast.util            - bounded nonlinear least squares, bivariate normal distribution
%
% Planned (released as the project progresses)
%   certaincoast.surrogate       - fast machine-learned copy of the model with a certified error bound
%   certaincoast.decisions       - safe operating window, pumping rule and safe pumping tests
%   (the model grows to compaction, liquefaction and coastal water levels in WP2)
%
% Start with run_tests and the scripts in the examples folder.
%
% Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.
