function st = lmStart(arm, rows, d, sigma, S, from)
%LMSTART Start a Levenberg-Marquardt calibration of one model structure (an arm) with a Gaussian prior.
%   ST = certaincoast.calibration.lmStart(ARM, ROWS, D, SIGMA) starts the calibration of the structure ARM on the
%   data D (the outputs ROWS of the model, independent noise of standard deviation SIGMA, a vector or a scalar) at
%   the prior mean ARM.th0, with one model run. ARM has the fields
%       fwd   handle: the column of all outputs of one model run, f = ARM.fwd(theta)
%       th0   prior mean of the parameters theta;   Gam   their prior covariance
%   The objective is the regularised least-squares (maximum a posteriori) objective
%       Phi(theta) = sum(((D - f(ROWS)) ./ SIGMA).^2) + (theta - th0)' Gam^-1 (theta - th0),
%   where the prior acts as Tikhonov regularisation, as in PEST.
%   ST = certaincoast.calibration.lmStart(ARM, ROWS, D, SIGMA, S, TH) starts at the parameters TH; S holds the
%   settings (certaincoast.calibration.defaults).
%   ST = certaincoast.calibration.lmStart(ARM, ROWS, D, SIGMA, S, ST0) continues the calibration ST0 on new data
%   from its current iterate with no new run: the outputs and the Jacobian of ST0 cover every output of the model.
%   It is used to calibrate the selected structure on all data at the end of the bandit.
%
%   ST has the fields th (iterate), f (all outputs at th), Phi, rows, d, w (weights 1/SIGMA.^2), th0, Gi (inverse
%   prior covariance), J and Jlast (finite-difference Jacobians at the iterate and at the previous iterate), lam
%   (damping), iter, conv, runs (model runs so far) and costLast (runs of the last call).
%
%   See also certaincoast.calibration.lmStep, certaincoast.calibration.lmPredict, certaincoast.calibration.trialGo.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 5 || isempty(S)
    S = certaincoast.calibration.defaults();
else
    S = certaincoast.calibration.defaults(S);
end
if ~isstruct(arm) || ~isscalar(arm) || ~all(isfield(arm, {'fwd', 'th0', 'Gam'}))
    error('certaincoast:calibration:badArm', 'ARM must be a structure with the fields fwd, th0 and Gam.');
end
rows = rows(:);
d = d(:);
sigma = sigma(:);
if isscalar(sigma)
    sigma = sigma + zeros(size(d));
end
if numel(rows) ~= numel(d) || numel(sigma) ~= numel(d) || isempty(d)
    error('certaincoast:calibration:badData', 'ROWS, D and SIGMA must have the same number of values.');
end
if any(~(sigma > 0))
    error('certaincoast:calibration:badSigma', 'SIGMA must be positive.');
end
th0 = arm.th0(:);
Gam = arm.Gam;
if ~isequal(size(Gam), [numel(th0) numel(th0)])
    error('certaincoast:calibration:badPrior', 'Gam must be a square matrix of the size of th0.');
end
Gi = inv(Gam);
st = struct('rows', rows, 'd', d, 'w', 1 ./ sigma .^ 2, 'th0', th0, 'Gi', 0.5 * (Gi + Gi.'));
if nargin >= 6 && isstruct(from)
    st.th = from.th;
    st.f = from.f;
    st.J = from.J;
    st.Jlast = from.Jlast;
    st.runs = 0;
else
    th = th0;
    if nargin >= 6 && ~isempty(from)
        th = from(:);
    end
    st.th = th;
    st.f = arm.fwd(th);
    st.f = st.f(:);
    st.J = [];
    st.Jlast = [];
    st.runs = 1;
end
if max(rows) > numel(st.f) || min(rows) < 1
    error('certaincoast:calibration:badRows', 'ROWS must index the outputs of the model (1..%d).', numel(st.f));
end
r = st.d - st.f(st.rows);
dt = st.th - st.th0;
st.Phi = sum(st.w .* r .^ 2) + dt.' * st.Gi * dt;
st.lam = S.lam0;
st.iter = 0;
st.conv = false;
st.costLast = st.runs;
end
