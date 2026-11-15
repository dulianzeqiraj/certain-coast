% certaincoast.calibration
%   Calibration when the model may be wrong: TRIAL-GO, a robust bandit over competing model structures.
%   Each structure (an arm) is calibrated by Levenberg-Marquardt with a prior, one iteration per play; every play is
%   checked against a validation item the structure has not been fitted to (for example a metered pumping stop) and
%   costs 1 when the item falls outside the structure's predictive interval. The robust Gittins index of Cohen and
%   Treetanthiploet (2022) decides which structure to play next, so the budget of model runs goes to the structures
%   that predict; the structure of smallest miss rate is selected and calibrated on all data.
%
% Bandit and robust index
%   trialGo            - TRIAL-GO over competing structures (also UCB, Thompson and greedy indices for comparison)
%   robustGittinsTable - robust Gittins index of an arm with costs 0 or 1 under ambiguity (Appendix D of the paper)
%   robustGittinsIndex - the index at a state (p, n), read from the table
%   defaults           - settings of the package
%
% Calibration of one structure
%   lmStart            - start (or continue on new data) a regularised Levenberg-Marquardt calibration
%   lmStep             - one iteration: finite-difference Jacobian and one trial run
%   lmPredict          - linearised (Laplace) predictive of any model output
%
% Synthetic coastal aquifer (certaincoast.model.planAquifer and planHeads)
%   aquiferOutputs     - heads under the current pumping, metered pumping stops and the goal, from one run
%   aquiferArms        - five competing structures: zones with recharge fixed, uniform, facies, zones, zones + strip
%   aquiferTruth       - a truth richer than every structure, and its data
%   aquiferGoal        - future heads at the supply wells and the risk that either is intruded
%
% Start with examples/example_trial_go.m.
%
% Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.
