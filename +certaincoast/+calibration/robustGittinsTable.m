function T = robustGittinsTable(k, beta, H, np, N, ng)
%ROBUSTGITTINSTABLE Robust Gittins index of an arm whose plays cost 1 (a miss) or 0, under ambiguity.
%   T = certaincoast.calibration.robustGittinsTable(K, BETA, H) tabulates the robust Gittins index of Cohen and
%   Treetanthiploet (2022, Section 5 and Appendix D) for an arm whose plays cost h = 1 (a miss) with an unknown
%   probability theta, and h = 0 otherwise. K in [0, 1) is the credible level of the set of miss rates the player
%   considers, BETA in (0, 1) the discount and H the number of plays looked ahead (after them the player stops).
%   T = certaincoast.calibration.robustGittinsTable(K, BETA, H, NP, N, NG) sets the number NP of nodes of the miss
%   rate p in [0, 1] (default 101), the counts n = 1..N (default N = 80; a vector N gives the counts wanted) and the
%   number NG of nodes of the index in [0, 1] (default 401).
%
%   State (p, n): p the posterior-mean miss rate, n the count (with p = (s + 1/2)/(n + 1) after n plays with s
%   misses, the state (p, n + 1) is the Jeffreys posterior Beta(s + 1/2, n - s + 1/2)). The set of miss rates is the
%   central credible interval of probability K of Beta(p n, (1 - p) n),
%       [p-, p+] = [B^-1(1/2 - K/2), B^-1(1/2 + K/2)],
%   and the one-step expectation is the worst case over it, E(f) = max over theta in {p-, p+} of
%   theta f(1) + (1 - theta) f(0) (a linear function of theta takes its maximum at an end). For a fee gamma per play,
%       V_H = 0,   V_t(p, n + t) = min{0, E_t[h - gamma + beta V_{t+1}(p_{t+1}, n + t + 1)]},
%       p_{t+1} = (p (n + t) + h)/(n + t + 1)   (linear interpolation in p),
%   the player may stop (value 0) at every play but the first, U(gamma, p, n) = E_0[h - gamma + beta V_1], and the
%   index is the smallest fee at which playing pays, gamma(p, n) = min{gamma : U(gamma, p, n) <= 0}, refined linearly
%   between the two nodes that bracket U = 0. The arm of smallest index is played. As n grows, gamma - p tends to 0
%   like 1/sqrt(n): an exploration bonus at small n (gamma < p), a premium for ambiguity that grows with K, and a
%   smaller index for a more patient player (larger BETA).
%
%   T has the fields k, beta, H, p (1 x NP), n (1 x numel(N)), gam (NP x numel(N)) and ng; read it with
%   certaincoast.calibration.robustGittinsIndex. The full default table takes a few seconds; it is computed once.
%
%   Example
%       T = certaincoast.calibration.robustGittinsTable(0.5, 0.95, 60, 101, [1 20 1000], 201);
%       interp1(T.p, T.gam, 0.3) - 0.3          % about -0.14, 0.02 and 0.005: tends to 0
%
%   Reference: S. N. Cohen and T. Treetanthiploet (2022), Gittins' theorem under uncertainty, Electronic Journal of
%   Probability 27, paper 17, https://doi.org/10.1214/22-EJP742
%
%   See also certaincoast.calibration.robustGittinsIndex, certaincoast.calibration.trialGo.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 4 || isempty(np), np = 101; end
if nargin < 5 || isempty(N), N = 80; end
if nargin < 6 || isempty(ng), ng = 401; end
if ~isscalar(k) || ~(k >= 0 && k < 1)
    error('certaincoast:calibration:badK', 'k must be a scalar in [0, 1).');
end
if ~isscalar(beta) || ~(beta > 0 && beta < 1)
    error('certaincoast:calibration:badBeta', 'beta must be a scalar in (0, 1).');
end
if ~isscalar(H) || H < 1 || H ~= round(H)
    error('certaincoast:calibration:badHorizon', 'H must be a positive integer.');
end
if np < 3 || ng < 3 || np ~= round(np) || ng ~= round(ng)
    error('certaincoast:calibration:badGrid', 'np and ng must be integers of at least 3.');
end
if isscalar(N)
    nlist = 1:N;
else
    nlist = N(:).';
end
if isempty(nlist) || any(nlist < 1) || any(nlist ~= round(nlist))
    error('certaincoast:calibration:badCounts', 'the counts n must be positive integers.');
end

p = linspace(0, 1, np);
G = linspace(0, 1, ng).';
ep = 1e-6;                                         % a shape p n or (1 - p) n equal to 0 is replaced by 1e-6
qlo = 0.5 - k / 2;
qhi = 0.5 + k / 2;
cmax = max(nlist) + H;
PL = zeros(cmax, np); PH = zeros(cmax, np);
IU = zeros(cmax, np); WU = zeros(cmax, np); ID = zeros(cmax, np); WD = zeros(cmax, np);
need = false(cmax, 1);
for n = nlist
    need(n:n + H - 1) = true;
end
cs = find(need);                                   % all counts needed, at once
a = max(cs * p, ep);
b = max(cs * (1 - p), ep);
PL(cs, :) = betaQuantile(qlo, a, b);
PH(cs, :) = betaQuantile(qhi, a, b);
[IU(cs, :), WU(cs, :)] = linearWeights(bsxfun(@rdivide, cs * p + 1, cs + 1), np);
[ID(cs, :), WD(cs, :)] = linearWeights(bsxfun(@rdivide, cs * p, cs + 1), np);
gam = zeros(np, numel(nlist));
for in = 1:numel(nlist)
    n = nlist(in);
    V = zeros(ng, np);                             % V_H = 0
    U = V;
    for t = H - 1:-1:0
        c = n + t;
        Vu = bsxfun(@times, V(:, IU(c, :)), WU(c, :)) + bsxfun(@times, V(:, IU(c, :) + 1), 1 - WU(c, :));
        Vd = bsxfun(@times, V(:, ID(c, :)), WD(c, :)) + bsxfun(@times, V(:, ID(c, :) + 1), 1 - WD(c, :));
        f1 = bsxfun(@plus, 1 - G, beta * Vu);
        f0 = bsxfun(@plus, -G, beta * Vd);
        E = max(bsxfun(@times, PL(c, :), f1) + bsxfun(@times, 1 - PL(c, :), f0), ...
                bsxfun(@times, PH(c, :), f1) + bsxfun(@times, 1 - PH(c, :), f0));
        if t > 0
            V = min(0, E);
        else
            U = E;
        end
    end
    for ip = 1:np                                  % U decreases in gamma
        j = find(U(:, ip) <= 0, 1);
        if isempty(j)
            gam(ip, in) = 1;
        elseif j == 1
            gam(ip, in) = G(1);
        else
            u0 = U(j - 1, ip);
            u1 = U(j, ip);
            gam(ip, in) = G(j - 1) + (G(j) - G(j - 1)) * u0 / (u0 - u1);
        end
    end
end
T = struct('k', k, 'beta', beta, 'H', H, 'p', p, 'n', nlist, 'gam', gam, 'ng', ng);
end

function q = betaQuantile(prob, a, b)
% quantile of Beta(a, b) by bisection on betainc (the same in MATLAB and GNU Octave, robust for tiny shapes)
if prob <= 0
    q = zeros(size(a));
    return
end
if prob >= 1
    q = ones(size(a));
    return
end
lo = zeros(size(a));
hi = ones(size(a));
for it = 1:42
    x = 0.5 * (lo + hi);
    up = betainc(x, a, b) < prob;
    lo(up) = x(up);
    hi(~up) = x(~up);
end
q = 0.5 * (lo + hi);
end

function [i, w] = linearWeights(x, np)
% node i and weight w of node i for linear interpolation at x in [0, 1] on np equal nodes
s = x * (np - 1);
i = floor(s) + 1;
i = min(max(i, 1), np - 1);
w = 1 - (s - (i - 1));
end
