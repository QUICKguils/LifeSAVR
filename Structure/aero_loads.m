% TODO:
% - Rework pos and d naming.
% - Rotation matrix with vpasolve.
% - Implement vpasolve more idiomatically.
% - Wait for more precise CAD values.

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
[rho, ~, ~, ~] = ISA(h);       % Air properties at desired altitude.
TAS2EAS = sqrt(rho/C.rho_sl);  % Conversion from true airspeed to equivalent airspeed.
CP = D.FE.CP;                  % Extract the CP table, just for conciseness.

% From the statement.
theta_dd = 1.0472;  % Additional pitch acceleration [rad/s²].
psi_max  = 15;      % Maximum yaw angle allowed [°].

% The return value consists of a table
% that contains the computed aerodynamic loads for all the CP.
Loads = table(...
	'Size', [height(CP), 6], ...
	'VariableNames', {'n',      'aoa',    'L',      'P',      'F_fin',  'M_fus'}, ...
	'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double'});

% Distances bewteen airplane COG and relevant components [m].
Plane_COG     = [D.Plane.COG_x, D.Plane.COG_z];
dist.wing2COG = D.Comp{"Wings",           ["COG_x", "COG_z"]} - Plane_COG;  % TODO: verify that MAC is at wing COG.
dist.HT2COG   = D.Comp{"Horizontal tail", ["COG_x", "COG_z"]} - Plane_COG;
dist.VT2COG   = D.Comp{"Vertical tail",   ["COG_x", "COG_z"]} - Plane_COG;
dist.T2COG    = D.Comp{"Engine",          ["COG_x", "COG_z"]} - Plane_COG;

% Solve aircraft dynamic equilibrium for each CP.
% Write the results in Loads.
for idx = 1:height(CP)
	n   = CP{idx, 'n'};
	EAS = CP{idx, 'EAS'};
	loads = solve_loads(n, EAS);  % TODO: check if syms defs and vpasolve are ok to loop.
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

		% Retrieve the true airspeed.
		TAS = EAS / TAS2EAS;

		% Thrust [N].
		T = D.Propu.T_sls * thrust_SLSconv(TAS, h, D.Propu.engine.BPR, D.Propu.engine.G);

		% Plane weight [N].
		W = D.Plane.MTOW * C.g;

		% Drag [N].  TODO: modify when we have more precise values.
		%
		% Drag on the wings. We take only the wing induced drag as a first
		% estimation.
		CD_wing = D.Wing.CL_cr^2 / (D.Wing.e * pi * D.Wing.AR);
		D_wing  = 0.5 * rho * TAS^2 * D.Wing.surf * CD_wing;
		% Drag on the fuselage. We take the parasitic drag at zero lift of the
		% entire plane, as a first estimation.
		CD_body = D.Plane.CD_0;
		D_body  = 0.5 * rho * TAS^2 * D.Wing.surf * CD_body;  % TODO: check if wing surf or plane wetted surf.

		% Pitching moment [N*m].
		% This formula does not apply to our geometry.
		% C_M_0 = -0.1;
		% K_n   = 0.09;
		% C_M   = C_M_0 - K_n * D.Wing.CL_max - 0.0015 * psi_max;
		cm = D.Wing.airfoil.cm;
		M  = 0.5 * rho * TAS^2 * D.Wing.surf * D.Wing.smc * cm;

		% Plane inertia along Y-axis [kg*m²].
		I_theta = 3000;  % TODO: to be computed with the CAD.

		% Angle of attack, wing lift and HT lift are coupled through three
		% equations:
		% - wing lift,
		% - vertical dynamic equilibrium, and
		% - moments equilibrium around the plane COG.
		% We solve this system numerically, thanks to vpasolve().
		aoa   = sym('aoa');
		L     = sym('L');
		P     = sym('P');
		dx_w  = sym('dx_w');
		dz_w  = sym('dz_w');
		dx_HT = sym('dx_HT');
		dz_b  = sym('dz_b');
		eqns = [ ...
			L == 0.5 * C.rho_cr * TAS^2 * D.Wing.surf * D.Wing.CL_alpha * aoa, ...
			L + P == n*W - T*sind(aoa - D.Wing.aoi), ...
			L*dx_w - T*dist.T2COG(2) + D_body*dz_b + D_wing*dz_w + P*dx_HT + M == I_theta * theta_dd, ...
			dx_w  == cosd(aoa) * dist.wing2COG(1) - sind(aoa) * dist.wing2COG(2), ...
			dz_w  == sind(aoa) * dist.wing2COG(1) + cosd(aoa) * dist.wing2COG(2), ...
			dx_HT == cosd(aoa) * dist.HT2COG(1)   - sind(aoa) * dist.HT2COG(2), ...
			dz_b  == 0, ...  % TODO: Verify if this is valid.
		];
		% Solve the system.
		res = vpasolve(eqns, [aoa, L, P, dx_w, dz_w, dx_HT, dz_b]);
		% Unpack the results of interest.
		aoa_res = double(res.aoa);
		L_res   = double(res.L);
		P_res   = double(res.P);

		% Lift curve slope of the vertical tail.
		VT_a = 5.5 * D.VT.AR / (D.VT.AR + 2);

		% Lift of the vertical tail.
		F_fin = 0.5 * rho * EAS^2 * VT_a * D.HT.surf * psi_max;

		% NOTE: no need to take into account the M_tail.

		% Fuselage bending moment.
		M_fus = F_fin * dist.VT2COG(1);

		% Return the computed loads.
		Sol = table(n, aoa_res, L_res, P_res, F_fin, M_fus);
	end
end