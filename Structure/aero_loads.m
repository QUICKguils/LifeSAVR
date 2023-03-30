% TODO:
% - Encode tail data in data.mat.
% - Wait for the CAD to correcct the component placements.

function Loads = aero_loads(h)
% LOADS_AERO  Aerodynamic loads.
%
% This function computes the aerodynamic loads exerting on the wings,
% fuselage and tail, for all the critical points of the flight envelope.
%
% This function implements the methodology and formulae that can be
% found at:
% Aircraft Structures>lesson 6>slides 11 to 25.
%
% Parameter:
%   h: double
%	  Altitude at which the flight envelope is desired, in meter.
% Argument:
%   AeroLoads: table
%	  Aerodynamic loads exerting on the wings, fuselage and tail, for
%	  all the critical points of the flight envelope.

%% Imports

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Constants and aircraft data.
C = load(fullfile(file_dir, "../constants.mat"));
D = load(fullfile(file_dir, "../data.mat"));

% Utilities.
addpath(genpath(fullfile(file_dir, "../Utils")));

%% Options setting

% Default altitude is the cruise altitude.
if ~nargin
	h  = C.h_cr;
end

%% Main solve

% Function global variables.
[rho, ~, ~, a] = ISA(h);       % Air properties at desired altitude.
TAS2EAS = sqrt(rho/C.rho_sl);  % Conversion from true airspeed to equivalent airspeed.
CP = D.FE.CP;                  % Extract the CP table, just for conciseness.

% From the statement.
theta_dd_add = 1.0472;  % Additional pitch acceleration [rad/s²].
psi_max = 15;           % Maximum yaw angle allowed [°].


% The return value consists of a table that contains the computed
% aerodynamic loads for all the CP.
Loads = table(...
	'Size', [height(CP), 6], ...
	'VariableNames', {'n',      'aoa',    'L',      'P',      'F_fin',  'M_fus'}, ...
	'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double'});

% Solve aircraft dynamic equilibrium for each CP.
% Write the results in Loads.
for idx = height(CP)
	n   = CP{idx, 'n'};
	EAS = CP{idx, 'EAS'};
	loads = solve_loads(n, EAS);
	Loads(idx, :) = loads;
end

%% Solve aircraft dynamic equilibrium for the given CP.

	function Sol = solve_loads(n, EAS)
		% SOLVE_LOADS  Solve aircraft dynamic equilibrium for the given CP.
		%
		% Arguments:
		%   n   (double) -- CP load factor.
		%   EAS (double) -- CP equivalent airspeed, in [m/s].
		% Return:
		%   Sol (1x6 table) -- Solution of the dynamic equilibrium.

		% Thrust [N].
		% WARN: check if TAS or EAS.
		% Conversion of equivalent airspeed (EAS) to true airspeed (TAS).
		M = EAS / TAS2EAS / a;
		T = thrust_sl2alt(M, h);

		% Weight [N].
		nW(idx) = D.Plane.MTOW * C.g * n;

		% Drag [N].
		% TODO: must depend on the envelope point.
		CD_wing = Wing.CD_0 + Wing.CL_cr^2 / (Wing.e * pi * Wing.AR);
		% C_D_wing = max([c_D_cruise_loiter, c_D_cruise_in,  c_D_cruise_eg]);
		D_w(idx) = 0.5 * CD_wing * EAS^2 * rho * D.Wing.S;
		% TODO: modify when we have more precise values.
		CD_body = 0.015;
		D_b(idx)  = 0.5 * CD_body * EAS^2 * rho * 2*pi*width_f^2/4;

		% Pitching moment.
		C_M_0 = -0.1;
		K_n   = 0.09; % TODO: Check if this value is OK.
		% This formula does not apply to our geometry.
		% C_M = C_M_0 - K_n * D.Wing.CL_max - 0.0015 * psi_max;
		C_M = D.Wing.airfoil.C_m;
		M(idx) = 0.5 * C_M * EAS^2 * rho * D.Wing.S * D.Wing.smc;

		I_theta = 3000; % A calculer avec le CAD
		c_l_alpha_plane = 6.38;
		alpha_L0_f = -0.0539;

		% Angle of attack, wing lift and HT load are coupled through three
		% equations. We solve them numerically, thanks to vpasolve().
		% TODO: change once we have the z_pos of the required components
		% (wait for the CAD).
		syms aoa L P;
		eqI   = c_l_alpha_plane * (aoa - alpha_L0_f) - L / (0.5 * C.rho_cr * D.Wing.S * (EAS^2)); %Lift of the wing
		eqII  = abs(c_g - pos(1)) * L - abs(z_cg - z_pos(8)) * T(idx) + z_cg * D_b(idx) - z_pos(1) * D_w(idx) - abs(pos(4) - c_g) * P + M(idx) - I_theta * theta_dd_add; % Moments about COG
		eqIII = L + P - nW(idx) + T(idx) * sin(aoa-alpha_L0_f); % Vertical equilibrium
		[aoa, L, P] = double(vpasolve(eqI, eqII, eqIII, aoa, L, P));

		% Lift curve slope of the vertical tail.
		a_VT = 5.5 * AR_VT / (AR_VT + 2);

		% Lift of the vertical tail.
		F_fin = 0.5 * rho * EAS^2 * a_VT* S_VT * psi_max;

		% NOTE: no need to take into account the M_tail.

		% Fuselage bending moment.
		M_fus = F_fin * (z_cg - z_pos(3));

		% Return the computed loads.
		Sol = table([n, aoa, L, P, F_fin, M_fus]);
	end
end