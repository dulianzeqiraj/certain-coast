function P = henryParameters()
%HENRYPARAMETERS Physical parameters of the Henry seawater-intrusion benchmark.
%   P = certaincoast.model.henryParameters() returns
%       eps   density ratio (rho_sea - rho_fresh) / rho_fresh            0.025
%       phi   porosity                                                    0.35
%       Dmol  pore-scale diffusivity (m2/s); the bulk coefficient in the
%             transport equation is phi * Dmol = 6.6e-6 m2/s              1.88571e-5
%       Q     freshwater inflow per unit height (m2/s)                    6.6e-5
%       Kref  reference hydraulic conductivity (m/s)                      1.0e-2
%   so that the Henry number is phi * Dmol / Q = 0.1 (Voss and Souza, 1987).

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

P = struct('eps', 0.025, 'phi', 0.35, 'Dmol', 1.62925 / 86400, 'Q', 6.6e-5, 'Kref', 1.0e-2);
end
