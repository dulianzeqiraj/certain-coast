function [gain, ess] = exactGoalEig(Yin, Gin, Yout, design, noiseSd, chunk)
%EXACTGOALEIG Exact information a design gives about yes/no goals, from an ensemble of model runs.
%   GAIN = certaincoast.network.exactGoalEig(YIN, GIN, YOUT, DESIGN, NOISESD) estimates, without any linearisation,
%   the mutual information between the data of the sensors DESIGN and each binary goal:
%       YIN    nIn-by-m   noise-free observations at the m candidate sensors for nIn prior realisations
%       GIN    nIn-by-q   the q goals (true/false) for the same realisations, e.g. "well A intruded"
%       YOUT   nOut-by-m  observations (with noise) of other, independent prior realisations: the possible data
%       NOISESD           standard deviation of the observation noise
%   For every possible data set the realisations are re-weighted by their likelihood, and GAIN(k) = H(P0_k) -
%   mean H(P_k | data), H the binary entropy (nats). [GAIN, ESS] also returns the mean effective sample size of the
%   weights: when it is small the inner ensemble is too small for that design.
%   [...] = certaincoast.network.exactGoalEig(..., CHUNK) processes CHUNK outer realisations at a time (default 250).
%
%   See also certaincoast.network.thresholdGoalEig.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 6
    chunk = 250;
end
G = double(Gin);
Yd = Yin(:, design);
Od = Yout(:, design);
nOut = size(Od, 1);
H = zeros(1, size(G, 2));
essSum = 0;
y2 = sum(Yd .^ 2, 2).';
for a = 1:chunk:nOut
    b = min(a + chunk - 1, nOut);
    O = Od(a:b, :);
    D2 = bsxfun(@plus, sum(O .^ 2, 2), y2) - 2 * O * Yd.';        % squared distances, outer-by-inner
    S = -D2 / (2 * noiseSd ^ 2);
    S = bsxfun(@minus, S, max(S, [], 2));
    w = exp(max(S, -80));
    sw = sum(w, 2);
    p = bsxfun(@rdivide, w * G, sw);
    H = H + sum(hbin(p), 1);
    essSum = essSum + sum(sw .^ 2 ./ sum(w .^ 2, 2));
end
gain = hbin(mean(G, 1)) - H / nOut;
ess = essSum / nOut;
end

function v = hbin(p)
p = min(max(p, 0), 1);
v = zeros(size(p));
k = p > 0 & p < 1;
v(k) = -(p(k) .* log(p(k)) + (1 - p(k)) .* log(1 - p(k)));
end
