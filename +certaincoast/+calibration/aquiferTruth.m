function tr = aquiferTruth(M, pump, seed, withStrip, S)
%AQUIFERTRUTH A synthetic truth of the aquifer, richer than every structure, and its data.
%   TR = certaincoast.calibration.aquiferTruth(M, PUMP, SEED, WITHSTRIP) draws, after rng(SEED), the true
%   log-conductivity multiplier field
%       log K - log K0c = theta_zone(zone) + zeta + stripLog [cell within the coastal strip, when WITHSTRIP],
%   with theta_zone ~ N(0, M.zonePrior) and zeta a sub-zonal Gaussian field (standard deviation S.zetaSd, covariance
%   exp(-distance / S.zetaL)) that no structure of aquiferArms can represent; the true recharge
%   M.recharge exp(xi_R), xi_R ~ N(0, S.sigR^2); the future change xi_f ~ N(0, S.sigXi^2) of the recharge factor;
%   and the data: heads at the observation wells under the current pumping (noise S.sigH) and the 2J metered head
%   changes of the stops at A and B (noise S.sigV).
%   TR = certaincoast.calibration.aquiferTruth(M, PUMP, SEED, WITHSTRIP, S) uses the settings S
%   (certaincoast.calibration.defaults).
%
%   TR has the fields seed, hasStrip, thz, zeta, xiR, xif, u, recharge, f (true outputs of aquiferOutputs),
%   dH (J heads), dV (2J metered changes), Hfut0 and hnat (future and natural heads at A and B), Hfut (future heads
%   with the realised xi_f) and riskTrue, the probability over xi_f that A or B is intruded (future head below M.hc).
%   Synthetic data only.
%
%   See also certaincoast.calibration.aquiferArms, certaincoast.calibration.aquiferGoal.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 5 || isempty(S)
    S = certaincoast.calibration.defaults();
else
    S = certaincoast.calibration.defaults(S);
end
rng(seed);
thz = chol(M.zonePrior, 'lower') * randn(M.Nz, 1);
w = randn(M.Nc, 1);
if S.zetaSd > 0
    zeta = subzonalFactor(M, S.zetaSd, S.zetaL) * w;
else
    zeta = zeros(M.Nc, 1);
end
xiR = S.sigR * randn;
xif = S.sigXi * randn;
J = numel(M.obsCell);
noiseH = S.sigH * randn(J, 1);
noiseV = S.sigV * randn(2 * J, 1);
tr = struct('seed', seed, 'hasStrip', logical(withStrip), 'thz', thz, 'zeta', zeta, 'xiR', xiR, 'xif', xif);
tr.u = thz(M.zone) + zeta + S.stripLog * double(logical(withStrip)) * double(M.strip);
tr.recharge = M.recharge * exp(xiR);
tr.f = certaincoast.calibration.aquiferOutputs(M, tr.u, tr.recharge, pump);
tr.dH = tr.f(1:J) + noiseH;
tr.dV = tr.f(J + 1:3 * J) + noiseV;
H0 = tr.f(3 * J + 1:3 * J + 2);
hn = tr.f(3 * J + 3:3 * J + 4);
tr.Hfut0 = H0;
tr.hnat = hn;
tr.Hfut = H0 + xif * hn;
tr.riskTrue = 0.5 * erfc(-(max((M.hc - H0) ./ hn) / S.sigXi) / sqrt(2));
end

function L = subzonalFactor(M, sd, ell)
% lower Cholesky factor of the covariance sd^2 exp(-d / ell) over the cells (cached for the session)
persistent Lc key
k = [M.Nc sd ell sum(M.X) sum(M.Y)];
if isempty(Lc) || ~isequal(key, k)
    D = sqrt(bsxfun(@minus, M.X, M.X.') .^ 2 + bsxfun(@minus, M.Y, M.Y.') .^ 2);
    C = sd ^ 2 * exp(-D / ell);
    C(1:M.Nc + 1:end) = C(1:M.Nc + 1:end) + 1e-9;
    Lc = chol(C, 'lower');
    key = k;
end
L = Lc;
end
