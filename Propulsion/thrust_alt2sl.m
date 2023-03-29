function [alt_over_sl] = thrust_alt2sl(M, p)
% thrust_sl2alt  Thrust conversion from sea level to h.

%% Imports

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Constants and aircraft data.
C = load(fullfile(file_dir, "../constants.mat"));
load(fullfile(file_dir, "../data.mat"), "Propu");
BPR = Propu.engine.BPR;
G = Propu.engine.G;

% Utilities.
addpath(genpath(fullfile(file_dir, "../Utils")));

% Preterms computation.
A = -0.4327 * (C.p_cr/C.p_sl)^2 + 1.3855 * C.p_cr/C.p_sl + 0.0472;
X =  0.1377 * (C.p_cr/C.p_sl)^2 - 0.4374 * C.p_cr/C.p_sl + 1.3003;
Z =  0.9106 * (C.p_cr/C.p_sl)^2 - 1.7736 * C.p_cr/C.p_sl + 1.8697;

alt_over_sl = A ...
	- Z * p/C.p_sl * (0.377*(1+BPR)*M) / sqrt((1+0.82*BPR)*G) ...
	+ X * p/C.p_sl * (0.23+0.19*sqrt(BPR)) * M^2;
end