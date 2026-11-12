function GamF = faciesPrior(rhoF, blocks, cross)
%FACIESPRIOR Prior covariance of a facies-stratified aquifer with one cross-facies correlation.
%   GAMF = certaincoast.network.faciesPrior(RHOF, BLOCKS) = D + RHOF * C, where D = blkdiag(BLOCKS{:}) holds the
%   within-facies covariances (one SPD block per facies, for example sands, gravels and clays) and C couples the
%   facies through their normalised block means (C = 1/sqrt(n_k n_l) between two cells of different facies k and l).
%   RHOF = 1 pools the facies (a facies-blind prior); RHOF = 0 makes them independent. GAMF is positive
%   semi-definite for every RHOF in [0, 1] when each block has a variance of at least 1 along its mean direction, as
%   correlation blocks (1 - r) I + r 11' do; otherwise, and for any coupling that breaks positivity, an error is
%   raised instead of returning an invalid covariance.
%   GAMF = certaincoast.network.faciesPrior(RHOF, BLOCKS, CROSS) uses a given cross-facies coupling CROSS (n-by-n,
%   symmetric, zero on the diagonal blocks).
%
%   See also certaincoast.network.poolingLoss.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

D = blkdiag(blocks{:});
n = size(D, 1);
if nargin < 3 || isempty(cross)
    sz = cellfun(@(b) size(b, 1), blocks);
    u = zeros(n, 1);
    mask = zeros(n);
    idx = 1;
    for k = 1:numel(blocks)
        r = idx:idx + sz(k) - 1;
        u(r) = 1 / sqrt(sz(k));
        mask(r, r) = 1;
        idx = idx + sz(k);
    end
    cross = (u * u') .* (1 - mask);
end
GamF = D + rhoF * cross;
GamF = 0.5 * (GamF + GamF');
if min(eig(GamF)) < -1e-10 * max(1, norm(GamF))
    error('certaincoast:network:notPositive', ...
          'The prior covariance is not positive semi-definite for rhoF = %g: rescale the blocks or the coupling.', rhoF);
end
end
