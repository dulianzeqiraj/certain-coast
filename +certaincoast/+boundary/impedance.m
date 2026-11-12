function Z = impedance(s, C, w, mu)
%IMPEDANCE Impedance of a boundary made of a Robin term and storage compartments.
%   Z = certaincoast.boundary.impedance(S, C, W, MU) = C + sum_m W(m) S / (S + MU(m)) for the complex frequencies S
%   (for example S = 1i * OMEGA, OMEGA in rad/day). Z is a column vector.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

s = s(:);
Z = C + sum(bsxfun(@times, w(:).', bsxfun(@rdivide, s, bsxfun(@plus, s, mu(:).'))), 2);
end
