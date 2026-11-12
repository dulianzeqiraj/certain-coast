function [kHat, kLo, kHi, openLo, openHi] = timeScaleInterval(times, x, kGrid, level)
%TIMESCALEINTERVAL Maximum-likelihood time scale of the indicator and its likelihood-ratio confidence interval.
%   [KHAT, KLO, KHI] = certaincoast.betweensamples.timeScaleInterval(TIMES, X, KGRID) returns the estimate of the rate
%   KAPPA (1/days) of a standardised record X sampled at TIMES, and the ends of its 95% likelihood-ratio interval,
%   searched over the range of the grid KGRID (increasing, positive).
%   [...] = certaincoast.betweensamples.timeScaleInterval(..., LEVEL) sets the confidence level.
%   [KHAT, KLO, KHI, OPENLO, OPENHI] also says whether an end of the interval is open, that is, the record cannot
%   exclude the end of the grid. With quarterly samples the upper end is usually open: the samples cannot tell a
%   chloride that changes in days from one that changes in weeks.
%
%   Only core functions are used (fminbnd, fzero, erfcinv).

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 4
    level = 0.95;
end
f = @(q) -certaincoast.betweensamples.ouLogLikelihood(exp(q), times, x);
q = log(kGrid(:));
ll = zeros(size(q));
for i = 1:numel(q)
    ll(i) = -f(q(i));
end
[~, i] = max(ll);
a = q(max(i - 1, 1));
b = q(min(i + 1, numel(q)));
[qh, fval] = fminbnd(f, a, b, optimset('TolX', 1e-6));
cut = -fval - 0.5 * (sqrt(2) * erfcinv(1 - level)) ^ 2;   % chi-square(1) quantile = (normal quantile)^2
g = @(qq) -f(qq) - cut;
openLo = g(q(1)) > 0;
openHi = g(q(end)) > 0;
if openLo
    kLo = kGrid(1);
else
    kLo = exp(fzero(g, [q(1), qh]));
end
if openHi
    kHi = kGrid(end);
else
    kHi = exp(fzero(g, [qh, q(end)]));
end
kHat = exp(qh);
end
