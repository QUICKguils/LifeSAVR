function [alt_over_sl] = thrust_alt2sl(TAS, h, BPR, G)
% thrust_sl2alt  Thrust conversion from sea level to h.

% Make Utils/ISA visible.
addpath(genpath(fullfile(fileparts(mfilename("fullpath")), "../Utils")));

% Compute air properties at altitude h and at sea level.
% NOTE: we don't load constants.mat in this function, for speed.
[~, p,    ~, a] = ISA(h);
[~, p_sl, ~, ~] = ISA(0);

% Retrieve Mach number, for the given flight condition (TAS, h).
M = TAS / a;

% Preterms computation.
A = -0.4327 * (p/C.p_sl)^2 + 1.3855 * p/p_sl + 0.0472;
X =  0.1377 * (p/C.p_sl)^2 - 0.4374 * p/p_sl + 1.3003;
Z =  0.9106 * (p/C.p_sl)^2 - 1.7736 * p/p_sl + 1.8697;

alt_over_sl = A ...
	- Z * p/C.p_sl * (0.377*(1+BPR)*M) / sqrt((1+0.82*BPR)*G) ...
	+ X * p/C.p_sl * (0.23+0.19*sqrt(BPR)) * M^2;
end