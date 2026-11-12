function [G, H] = adjointSensitivities(J, dRdtheta, dOdu, dqdu)
%ADJOINTSENSITIVITIES Sensitivities of all candidate observations and of the goal from two adjoint solves.
%   [G, H] = certaincoast.network.adjointSensitivities(J, DRDTHETA, DODU, DQDU) for the Jacobian J = dR/du of the
%   model at its solution, the parameter derivative DRDTHETA = dR/dtheta, the observation operator DODU (one row per
%   candidate sensor) and the goal derivative DQDU (one row per goal) solves
%
%       J' Lambda = DODU',   J' mu = DQDU'     (one factorisation of J', all right-hand sides together)
%
%   and returns G = -Lambda' DRDTHETA (candidates-by-parameters) and H = -mu' DRDTHETA (goals-by-parameters).
%   The cost does not grow with the number of parameters, so a whole field of log-conductivities is affordable.
%
%   See also certaincoast.model.densityJacobian, certaincoast.network.eigOperator.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

m = size(dOdu, 1);
Z = J' \ [dOdu', dqdu'];
G = -(Z(:, 1:m)' * dRdtheta);
H = -(Z(:, m + 1:end)' * dRdtheta);
if issparse(G)
    G = full(G);
    H = full(H);
end
end
