function I = pairInformation(kappa, tau)
%PAIRINFORMATION Fisher information on the time scale from one pair of samples TAU days apart.
%   I = certaincoast.betweensamples.pairInformation(KAPPA, TAU) is the Fisher information on KAPPA carried by one
%   transition of a standardised stationary Ornstein-Uhlenbeck indicator over the lag TAU (Proposition S2 of the
%   preparatory experiment PATH-GO):
%
%       I = TAU^2 RHO^2 (1 + RHO^2) / (1 - RHO^2)^2,   RHO = exp(-KAPPA TAU).
%
%   It tends to 1 / (2 KAPPA^2) for short lags and decays as 2 TAU^2 exp(-2 KAPPA TAU) for long ones, so quarterly
%   samples say almost nothing about a chloride that changes in days. With n transitions, the relative standard error
%   of KAPPA cannot fall below 1 / (KAPPA sqrt(n I)).

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

r = exp(-kappa .* tau);
I = tau .^ 2 .* r .^ 2 .* (1 + r .^ 2) ./ (1 - r .^ 2) .^ 2;
end
