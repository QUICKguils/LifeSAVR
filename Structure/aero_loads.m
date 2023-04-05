% TODO:
% - Check aoi of the HT tail, it's probably different that for the wing.
% - D.Propu.T_sls should not be used for every CP.
% - Rework distances naming convention.
% - Try to define syms variables in the main function scope, for speed.
% - Implement vpasolve more idiomatically. Maybe rotation matrix for distances.
% - Wait for more precise CAD values.

function aero_loads
% LOADS_AERO  Aerodynamic loads.
%
% This function computes the aerodynamic loads exerting on the wings,
% fuselage and tail, for all the critical points of the flight envelope.
%
% This function implements the methodology and formulae that can be
% found at:
% Aircraft Structures>lesson 6>slides 11 to 25.
%
% Save:
%   AeroLoads: 5x7 table
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

%% Main solve

% Function global variables.
h = D.FE.altitude;             % Atitude at which the CP are computed [m].
[rho, ~, ~, ~] = ISA(h);       % Air properties at desired altitude.
TAS2EAS = sqrt(rho/C.rho_sl);  % Conversion from true airspeed to equivalent airspeed.
CP = D.FE.CP;                  % Extract the CP table, just for conciseness.

% From the statement.
theta_dd = 1.0472;       % Additional pitch acceleration [rad/s²].
psi_max  = deg2rad(15);  % Maximum yaw angle allowed [rad].

% The data structure to save consists of a table
% that contains the computed aerodynamic loads for all the CP.
AeroLoads = table(...
	'Size', [height(CP), 7], ...
	'VariableNames', {'n',      'EAS',    'aoa',    'L',      'P',      'F_fin',  'M_fus'}, ...
	'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double', 'double'});

% Distances bewteen airplane COG and relevant components [m].
dist.wing2COG = D.Comp{"Wings",           "COG"} - D.Plane.COG;  % TODO: verify that MAC is at wing COG.
dist.HT2COG   = D.Comp{"Horizontal tail", "COG"} - D.Plane.COG;
dist.E2COG    = D.Comp{"Engine",          "COG"} - D.Plane.COG;

% Solve aircraft dynamic equilibrium for each CP.
% Write the results in AeroLoads.
for idx = 1:height(CP)
	n   = CP{idx, 'n'};
	EAS = CP{idx, 'EAS'};
	loads = solve_loads(n, EAS);
	AeroLoads(idx, :) = loads;
end

% Save AeroLoads in data.mat.
save(fullfile(file_dir, "../data.mat"), "AeroLoads", "-append");

%% Solve aircraft dynamic equilibrium for the given CP

	function loads = solve_loads(n, EAS)
		% SOLVE_LOADS  Solve aircraft dynamic equilibrium for the given CP.
		%
		% Arguments:
		%   n   (double) -- CP load factor.
		%   EAS (double) -- CP equivalent airspeed, in [m/s].
		% Return:
		%   loads (1x6 table) -- Solution of the dynamic equilibrium.

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
		% Only take into account the wing pitching moment.
		cm = D.Wing.airfoil.cm;
		M  = 0.5 * rho * TAS^2 * D.Wing.surf * D.Wing.smc * cm;

		% Plane inertia along Y-axis [kg*m²].
		I_theta = sum(D.Comp.Mass .* vecnorm(D.Comp.COG - D.Plane.COG, 2, 2).^2);

		% Angle of attack, wing lift and HT lift are coupled through three
		% equations:
		% - wing lift,
		% - vertical dynamic equilibrium, and
		% - moments equilibrium around the plane COG.
		% Note that, for the moments equilibrium, the lever arms depend on
		% the angle of attack. This adds 4 more equations to the system.
		%
		% Define symbolic variables.
		aoa   = sym('aoa');    % Angle of attack [°].
		aoi   = sym('aoi');    % Angle of incidence [°].
		L     = sym('L');      % Lift from the wing [N].
		P     = sym('P');      % Lift from the horizontal tail [N].
		dx_w  = sym('dx_w');   % X-distance from airplaine COG to wing COG, absolute axis [m].
		dz_w  = sym('dz_w');   % Z-distance from airplaine COG to wing COG, absolute axis [m].
		dx_HT = sym('dx_HT');  % X-distance from airplaine COG to HT   COG, absolute axis [m].
		dz_b  = sym('dz_b');   % Z-distance from airplaine COG to body COG, absolute axis [m].
		% Define system of equations.
		eqns = [ ...
			L == 0.5 * rho * TAS^2 * D.Wing.surf * D.Wing.CL_alpha * deg2rad(aoa + D.Wing.aoa_zerolift), ...
			L + P == n*W - T*sind(aoi), ...
			I_theta * theta_dd == - L*dx_w + D_wing*dz_w + D_body*dz_b - P*dx_HT - T*dist.E2COG(3) + M, ...
			dx_w  ==   cosd(aoi) * dist.wing2COG(1) + sind(aoi) * dist.wing2COG(3), ...
			dz_w  == - sind(aoi) * dist.wing2COG(1) + cosd(aoi) * dist.wing2COG(3), ...
			dx_HT ==   cosd(aoi) * dist.HT2COG(1)   + sind(aoi) * dist.HT2COG(3), ...
			dz_b  == 0, ...  % TODO: Verify if this is valid.
			aoi == aoa - D.Wing.aoi, ...
		];
		% Solve the system.
		res = vpasolve(...
			eqns, ...
			[aoa, L,   P, dx_w,             dz_w,             dx_HT,          dz_b, aoi], ...
			[5,   n*W, 0, dist.wing2COG(1), dist.wing2COG(3), dist.HT2COG(1), 0,    5]);
		% Unpack the results of interest.
		aoa_res = double(res.aoa);
		L_res   = double(res.L);
		P_res   = double(res.P);

		% Lift curve slope of the vertical tail [1/rad].  % TODO: make sure of units.
		VT_a = 5.5 * D.VT.AR / (D.VT.AR + 2);

		% Lift of the vertical tail [N].
		F_fin = 0.5 * rho * EAS^2 * D.VT.surf * VT_a * psi_max;

		% Fuselage bending moment [N*m].
		% No need to take into account the M_tail.  % TODO: not sure of that: M_tail does not seems negligible.
		M_fus = F_fin * D.VT.y;  % TODO: verify that using D.VT.y is correct.

		% Return the computed loads.
		loads = table(n, EAS, aoa_res, L_res, P_res, F_fin, M_fus);
	end
end