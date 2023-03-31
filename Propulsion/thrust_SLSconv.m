function [SLSconv] = thrust_SLSconv(TAS, h, BPR, G)
% THRUST_SLSCONV  Conversion factor from sea level, static (SLS) thrust.
%
% This function provides the factor to make the conversion between the
% static thrust at sea level, and the equivalent thrust at the given
% flight conditions (true airspeed and altitude).
% 
% This function simply implements the conversion formulae given in:
% Teams>Propulsion>Files>Sources>Bartel-SFC_Calculations.pdf (page 8)
%
% Arguments:
%   TAS (double) -- True airspeed of the plane [m/s].
%   h   (double) -- Flight altitude [m].
%   BPR (double) -- Engine bypass ratio.
%   G   (double) -- Gas generator function.
% Return:
%   alt_over_sl (double)  -- Thrust conversion factor.

% Make Utils/ISA visible.
addpath(genpath(fullfile(fileparts(mfilename("fullpath")), "../Utils")));

% NOTE:
% This function don't load constants.mat. We keep it independent of any
% global data, for efficiency reasons. Indeed, loading MAT files is
% too penalizing for such a small routine.

% Compute air properties at altitude h and at sea level.
[~, p,    ~, a] = ISA(h);
[~, p_sl, ~, ~] = ISA(0);

% Retrieve Mach number, for the given flight condition (TAS, h).
M = TAS / a;

% Preterms computation.
A = -0.4327 * (p/p_sl)^2 + 1.3855 * p/p_sl + 0.0472;
X =  0.1377 * (p/p_sl)^2 - 0.4374 * p/p_sl + 1.3003;
Z =  0.9106 * (p/p_sl)^2 - 1.7736 * p/p_sl + 1.8697;

SLSconv = A ...
	- Z * p/p_sl * (0.377*(1+BPR)*M) / sqrt((1+0.82*BPR)*G) ...
	+ X * p/p_sl * (0.23+0.19*sqrt(BPR)) * M^2;
end