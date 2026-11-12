function [t, r, Q] = forcing(p, nDays, kind, drought, newPumping)
%FORCING Daily recharge and pumping for the test aquifer.
%   [T, R, Q] = certaincoast.boundary.forcing(P, NDAYS) returns times (days), recharge R (m/d, seasonal around P.R0)
%   and pumping Q (m2/d per metre of coast, constant P.Q0).
%   [...] = certaincoast.boundary.forcing(P, NDAYS, 'broadband') makes the pumping change at random every 6 months to
%   3 years, with a monthly modulation: a richer forcing that lets two piezometers identify the boundary.
%   [...] = certaincoast.boundary.forcing(P, NDAYS, KIND, DROUGHT, NEWPUMPING) multiplies the recharge by
%   DROUGHT(3) for DROUGHT(1) <= t < DROUGHT(2), and sets the pumping to NEWPUMPING(3) for NEWPUMPING(1) <= t <
%   NEWPUMPING(2). Use [] to skip either.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

if nargin < 3 || isempty(kind)
    kind = 'seasonal';
end
t = (1:nDays).' * p.dt;
r = p.R0 * (1 + 0.8 * sin(2 * pi * t / 365.25));
if nargin >= 4 && ~isempty(drought)
    in = t >= drought(1) & t < drought(2);
    r(in) = drought(3) * r(in);
end
switch kind
    case 'seasonal'
        Q = p.Q0 * ones(nDays, 1);
    case 'broadband'
        Q = zeros(nDays, 1);
        k = 0;
        while k < nDays
            n = 180 + floor(901 * rand);
            Q(k + 1:min(k + n, nDays)) = 0.3 + 0.7 * rand;
            k = k + n;
        end
        month = floor(t / 30.4375) + 1;
        f = 1 + 0.15 * (2 * rand(max(month), 1) - 1);
        Q = Q .* f(month);
    otherwise
        error('certaincoast:boundary:badKind', 'kind must be ''seasonal'' or ''broadband''.');
end
if nargin >= 5 && ~isempty(newPumping)
    Q(t >= newPumping(1) & t < newPumping(2)) = newPumping(3);
end
end
