function ll = ouLogLikelihood(kappa, times, x)
%OULOGLIKELIHOOD Exact log-likelihood of a standardised stationary Ornstein-Uhlenbeck record at irregular times.
%   LL = certaincoast.betweensamples.ouLogLikelihood(KAPPA, TIMES, X) for samples X (standardised: mean 0,
%   variance 1) taken at TIMES (days).

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

times = times(:);
x = x(:);
dtau = diff(times);
r = exp(-kappa * dtau);
v = 1 - r .^ 2;
res = x(2:end) - r .* x(1:end - 1);
ll = -0.5 * log(2 * pi) - 0.5 * x(1) ^ 2 - 0.5 * sum(log(2 * pi * v) + res .^ 2 ./ v);
end
