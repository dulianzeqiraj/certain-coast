function out = soilWaterBalance(par, Pg, PET)
%SOILWATERBALANCE Daily chain from gauge rain to recharge and to the head of one well.
%   OUT = certaincoast.recharge.soilWaterBalance(PAR, PG, PET) runs, day by day,
%
%       gauge rain PG --(undercatch m)--> rain P --(infiltration capacity Ic)--> infiltration
%       --(root-zone store theta, field capacity theta_fc, evapotranspiration k_et * PET * min(1, theta/theta_fc))-->
%       drainage (theta - theta_fc)+ / tau_s = recharge R --(specific yield Sy, drainage time tau_a)--> head h.
%
%   PAR is a struct with fields m, k_et, Ic (mm/d), theta_fc (mm), tau_s (d), phi0 (initial store as a fraction of
%   theta_fc), Sy, tau_a (d), h_b (m, drainage base) and h_0 (m, initial head). PG and PET are in mm/d.
%   OUT has fields h (m), R (mm/d), theta (mm), ET (mm/d) and runoff (mm/d), all column vectors.
%
%   The heads do not change when m, k_et, Ic, theta_fc and Sy are multiplied by the same factor (Proposition 2 of the
%   preparatory experiment INFILTRATION-GO): heads and rain alone cannot fix the amount of recharge; one anchor
%   (a pumping test for Sy, a soil-moisture probe or a chloride mass balance) is needed.
%
%   See also certaincoast.recharge.rechargeFromHeads.

%   Copyright 2026 Dulian Zeqiraj. Licensed under the Apache License, Version 2.0.

N = numel(Pg);
th = par.phi0 * par.theta_fc;
h = par.h_0;
out = struct('h', zeros(N, 1), 'R', zeros(N, 1), 'theta', zeros(N, 1), 'ET', zeros(N, 1), 'runoff', zeros(N, 1));
for i = 1:N
    P = par.m * Pg(i);
    inf_ = min(P, par.Ic);
    th = th + inf_;
    et = min(par.k_et * PET(i) * min(1, th / par.theta_fc), th);
    th = th - et;
    d = max(th - par.theta_fc, 0) / par.tau_s;
    th = th - d;
    h = h + d / (1000 * par.Sy) - (h - par.h_b) / par.tau_a;
    out.h(i) = h;
    out.R(i) = d;
    out.theta(i) = th;
    out.ET(i) = et;
    out.runoff(i) = P - inf_;
end
end
