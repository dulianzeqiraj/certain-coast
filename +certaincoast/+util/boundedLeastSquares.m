function [x, info] = boundedLeastSquares(fun, x0, lb, opts)
%BOUNDEDLEASTSQUARES Nonlinear least squares with lower bounds (projected Levenberg-Marquardt).
%   X = certaincoast.util.boundedLeastSquares(FUN, X0, LB) minimises 0.5 * sum(FUN(X).^2) subject to X >= LB, starting
%   from X0. FUN must accept a matrix whose columns are parameter vectors and return a matrix whose columns are the
%   residual vectors: the forward-difference Jacobian is then computed in a single call.
%   X = certaincoast.util.boundedLeastSquares(..., OPTS) with fields MaxIter (default 200), TolFun (relative decrease
%   of the cost, default 1e-10) and Step (relative difference step, default 1e-6).
%   [X, INFO] also returns the final cost, the number of iterations and the exit reason.
%
%   Only core functions are used, so no toolbox is needed, and the same code runs in MATLAB and GNU Octave.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 4
    opts = struct();
end
maxIter = getopt(opts, 'MaxIter', 200);
tolFun = getopt(opts, 'TolFun', 1e-10);
step = getopt(opts, 'Step', 1e-6);
lb = lb(:);
x = max(x0(:), lb);
n = numel(x);
r = fun(x);
cost = 0.5 * (r.' * r);
lambda = 1e-3;
reason = 'maximum number of iterations';
for it = 1:maxIter
    hstep = step * max(abs(x), 1e-3 * max(abs(x0(:)), 1));
    Xp = repmat(x, 1, n) + diag(hstep);
    J = bsxfun(@rdivide, bsxfun(@minus, fun(Xp), r), hstep.');
    g = J.' * r;
    A = J.' * J;
    free = ~(x <= lb & g > 0);                       % variables held at their bound
    improved = false;
    while lambda < 1e12
        dx = zeros(n, 1);
        Af = A(free, free);
        D = diag(max(diag(Af), eps));
        dx(free) = -(Af + lambda * D) \ g(free);
        xn = max(x + dx, lb);
        rn = fun(xn);
        cn = 0.5 * (rn.' * rn);
        if cn < cost
            improved = true;
            rel = (cost - cn) / max(cost, realmin);
            x = xn;
            r = rn;
            cost = cn;
            lambda = max(lambda / 3, 1e-12);
            break
        end
        lambda = lambda * 4;
    end
    if ~improved
        reason = 'no further decrease';
        break
    end
    if rel < tolFun
        reason = 'relative decrease below TolFun';
        break
    end
end
info = struct('cost', cost, 'iterations', it, 'reason', reason);
end

function v = getopt(opts, name, default)
if isfield(opts, name)
    v = opts.(name);
else
    v = default;
end
end
