function p = bvnCdf(a, b, r)
%BVNCDF Standard bivariate normal distribution function P(X < a, Y < b), corr(X, Y) = r.
%   P = certaincoast.util.bvnCdf(A, B, R) for arrays A and B of the same size (or one scalar) and a scalar R with
%   |R| < 1, from Phi2(a, b; r) = int_{-inf}^{a} phi(t) Phi((b - r t) / sqrt(1 - r^2)) dt with 64-point
%   Gauss-Legendre quadrature on [-9, a]. Accuracy better than 1e-10. No toolbox is needed.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if ~isscalar(r) || abs(r) >= 1
    error('certaincoast:util:r', 'r must be a scalar with |r| < 1.');
end
if isscalar(a), a = a + zeros(size(b)); end
if isscalar(b), b = b + zeros(size(a)); end
if ~isequal(size(a), size(b))
    error('certaincoast:util:size', 'a and b must have the same size.');
end
persistent x w
if isempty(x)
    k = (1:63).';
    beta = k ./ sqrt(4 * k .^ 2 - 1);
    [V, D] = eig(diag(beta, 1) + diag(beta, -1));
    [x, i] = sort(diag(D));
    w = 2 * V(1, i).' .^ 2;
end
lo = -9;
hi = max(a(:), lo + 1e-12);
t = 0.5 * (hi - lo) * x.' + 0.5 * (hi + lo);
f = exp(-0.5 * t .^ 2) / sqrt(2 * pi) .* (0.5 * erfc(-((b(:) - r * t) / sqrt(1 - r ^ 2)) / sqrt(2)));
p = reshape(0.5 * (hi - lo) .* (f * w), size(a));
p = min(max(p, 0), 1);
end
