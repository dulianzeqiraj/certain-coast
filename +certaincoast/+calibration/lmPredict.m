function [mu, C, Sig] = lmPredict(st, rows)
%LMPREDICT Linearised (Laplace) predictive of model outputs at the current iterate of a calibration.
%   [MU, C] = certaincoast.calibration.lmPredict(ST, ROWS) returns the mean MU = f(ROWS) at the current iterate of
%   the calibration ST (lmStart, lmStep) and the covariance C = J_r Sig J_r' of the outputs ROWS, with the
%   posterior covariance of the parameters
%       Sig = (J_d' W J_d + Gam^-1)^-1,
%   J the most recent finite-difference Jacobian (at the iterate if it is available, otherwise at the previous
%   iterate, as PEST does), J_d its rows of the data and J_r its rows ROWS. The noise of a new reading is not
%   included (add its variance to the diagonal of C), nor any error of the structure itself: a wrong structure is
%   detected by the misses of its intervals (certaincoast.calibration.trialGo).
%   [MU, C, SIG] = ... also returns Sig.
%
%   See also certaincoast.calibration.lmStart, certaincoast.calibration.lmStep.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

J = st.J;
if isempty(J)
    J = st.Jlast;
end
if isempty(J)
    error('certaincoast:calibration:noJacobian', 'no Jacobian yet: call lmStep once first.');
end
rows = rows(:);
Jd = J(st.rows, :);
Hm = Jd.' * bsxfun(@times, st.w, Jd) + st.Gi;
Hm = 0.5 * (Hm + Hm.');
Sig = inv(Hm);
Sig = 0.5 * (Sig + Sig.');
mu = st.f(rows);
Jr = J(rows, :);
C = Jr * Sig * Jr.';
C = 0.5 * (C + C.');
end
