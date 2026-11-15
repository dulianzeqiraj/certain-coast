function st = lmStep(arm, st, S)
%LMSTEP One Levenberg-Marquardt iteration of a calibration started by lmStart (as one iteration of PEST).
%   ST = certaincoast.calibration.lmStep(ARM, ST) computes the forward-difference Jacobian at the current iterate
%   (n model runs for n parameters, step S.fdStep; kept when it is still valid after a rejected step) and makes one
%   trial run at
%       theta + dtheta,   (J' W J + Gam^-1 + lam diag(J' W J + Gam^-1)) dtheta = J' W r - Gam^-1 (theta - th0),
%   with r the residuals of the data and W = diag(1/sigma^2). The step is accepted when it lowers the objective
%   (the damping lam is then divided by 4), otherwise lam is multiplied by 8. The calibration has converged when an
%   accepted step lowers the objective by less than the fraction S.tolRel, when lam exceeds 1e6 or after S.maxIt
%   iterations; a call after convergence costs no run. ST.costLast holds the runs of this call.
%   ST = certaincoast.calibration.lmStep(ARM, ST, S) uses the settings S (certaincoast.calibration.defaults).
%
%   See also certaincoast.calibration.lmStart, certaincoast.calibration.lmPredict.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 3 || isempty(S)
    S = certaincoast.calibration.defaults();
else
    S = certaincoast.calibration.defaults(S);
end
st.costLast = 0;
if st.conv
    return
end
n = numel(st.th);
c = 0;
if isempty(st.J)
    J = zeros(numel(st.f), n);
    for i = 1:n
        tp = st.th;
        tp(i) = tp(i) + S.fdStep;
        fp = arm.fwd(tp);
        J(:, i) = (fp(:) - st.f) / S.fdStep;
    end
    c = c + n;
    st.J = J;
end
Jd = st.J(st.rows, :);
r = st.d - st.f(st.rows);
Hm = Jd.' * bsxfun(@times, st.w, Jd) + st.Gi;
Hm = 0.5 * (Hm + Hm.');
g = Jd.' * (st.w .* r) - st.Gi * (st.th - st.th0);
dth = (Hm + st.lam * diag(diag(Hm))) \ g;
tt = st.th + dth;
ft = arm.fwd(tt);
ft = ft(:);
c = c + 1;
rt = st.d - ft(st.rows);
dt = tt - st.th0;
Pt = sum(st.w .* rt .^ 2) + dt.' * st.Gi * dt;
if Pt < st.Phi
    rel = (st.Phi - Pt) / st.Phi;
    st.Jlast = st.J;
    st.J = [];
    st.th = tt;
    st.f = ft;
    st.Phi = Pt;
    st.lam = max(st.lam / 4, 1e-7);
    if rel < S.tolRel
        st.conv = true;
    end
else
    st.lam = st.lam * 8;
    if st.lam > 1e6
        st.conv = true;
    end
end
st.iter = st.iter + 1;
if st.iter >= S.maxIt
    st.conv = true;
end
st.runs = st.runs + c;
st.costLast = c;
end
