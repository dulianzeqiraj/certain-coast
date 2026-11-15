function H = planHeads(M, u, recharge, Q)
%PLANHEADS Steady heads of the plan-view confined aquifer of planAquifer, from one factorisation.
%   H = certaincoast.model.planHeads(M, U, RECHARGE) returns the Nc-by-3 matrix [hR, gA, gB]: hR is the head [m]
%   under the uniform recharge RECHARGE [m/d] alone, gA and gB the head change [m per m3/d] due to the extraction
%   of 1 m3/d at the supply wells A and B (negative: pumping lowers the head). U is the Nc-by-1 field of natural-log
%   multipliers of the conductivity, K = M.K0c .* exp(U); U = 0 is the background field. The head under any pumping
%   of the supply wells is the combination hR + QA gA + QB gB.
%   H = certaincoast.model.planHeads(M, U, RECHARGE, Q), with Q an Nc-by-m matrix of extra extraction rates [m3/d]
%   (injection: negative), appends the m columns of the head change due to each column of Q.
%
%   Finite volumes on square cells: the conductance of a face is the harmonic mean of the transmissivities K b of
%   its two cells; the coastal cells hold the sea level (head 0) and the other edges are closed. The matrix of the
%   other cells is symmetric positive definite and is factorised once (sparse Cholesky), so every pumping state is
%   one more right-hand side of the same run. Two consequences are used by certaincoast.calibration:
%     reciprocity  the head change at a well j per unit extraction at A equals the head change at A per unit
%                  extraction at j (the matrix is symmetric), so a metered pumping stop at A measures the
%                  sensitivity of the head at A to every well that is read;
%     invariance   multiplying every conductivity and the recharge by the same factor c leaves hR unchanged and
%                  divides gA and gB by c: heads under recharge cannot fix the scale of the fluxes, a metered
%                  change of pumping can.
%
%   Example
%       M = certaincoast.model.planAquifer();
%       H = certaincoast.model.planHeads(M, zeros(M.Nc, 1), M.recharge);
%       sA = -5000 * H(M.supplyCell(1), 2);      % drawdown at A for 5000 m3/d extracted at A
%
%   See also certaincoast.model.planAquifer, certaincoast.calibration.aquiferOutputs.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 4 || isempty(Q)
    Q = zeros(M.Nc, 0);
end
u = full(u(:));
if numel(u) ~= M.Nc
    error('certaincoast:model:badField', 'U must have one value per cell (%d).', M.Nc);
end
if size(Q, 1) ~= M.Nc
    error('certaincoast:model:badPumping', 'Q must have one row per cell (%d).', M.Nc);
end
if ~isscalar(recharge) || ~isfinite(recharge)
    error('certaincoast:model:badRecharge', 'RECHARGE must be a finite scalar (m/d).');
end
T = M.K0c .* exp(u) * M.b;                                   % transmissivity, m2/d
Tf = 2 * T(M.ek) .* T(M.el) ./ (T(M.ek) + T(M.el));         % conductance of every face (square cells)
d = accumarray([M.ek; M.el], [Tf; Tf], [M.Nc 1]);
ii = M.interior;
nI = numel(ii);
p = M.pos;
f = M.faceInterior;
A = sparse([p(M.ek(f)); p(M.el(f)); (1:nI).'], [p(M.el(f)); p(M.ek(f)); (1:nI).'], ...
           [-Tf(f); -Tf(f); d(ii)], nI, nI);
B = zeros(nI, 3 + size(Q, 2));
B(:, 1) = recharge * M.cellArea;
B(p(M.supplyCell(1)), 2) = -1;
B(p(M.supplyCell(2)), 3) = -1;
B(:, 4:end) = -full(Q(ii, :));
[L, flag, perm] = chol(A, 'lower', 'vector');
if flag ~= 0
    error('certaincoast:model:notPositive', 'the flow matrix is not positive definite.');
end
X = zeros(size(B));
X(perm, :) = L.' \ (L \ B(perm, :));
H = zeros(M.Nc, size(B, 2));
H(ii, :) = X;
end
