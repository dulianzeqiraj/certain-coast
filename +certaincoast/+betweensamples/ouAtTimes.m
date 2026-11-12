function x = ouAtTimes(times, kappa)
%OUATTIMES A standardised stationary Ornstein-Uhlenbeck record observed at irregular times.
%   X = certaincoast.betweensamples.ouAtTimes(TIMES, KAPPA) simulates the record (column vector) at the given times.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

times = times(:);
x = zeros(numel(times), 1);
x(1) = randn;
for i = 2:numel(times)
    r = exp(-kappa * (times(i) - times(i - 1)));
    x(i) = r * x(i - 1) + sqrt(1 - r ^ 2) * randn;
end
end
