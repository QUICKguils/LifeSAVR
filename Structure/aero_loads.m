% TODO:
% - Encode tail data in data.mat.
% - Use comp position stored in data.mat.
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
[rho, ~, ~, ~] = ISA(h);       % Air properties at desired altitude.
TAS2EAS = sqrt(rho/C.rho_sl);  % Conversion from true airspeed to equivalent airspeed.
CP = D.FE.CP;                  % Extract the CP table, just for conciseness.

% From the statement.
theta_dd_add = 1.0472;  % Additional pitch acceleration [rad/s²].
psi_max = 15;           % Maximum yaw angle allowed [°].


% The return value consists of a table
% that contains the computed aerodynamic loads for all the CP.
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

		% Retrieve the true airspeed.
		TAS = EAS / TAS2EAS;

		% Thrust [N].
		T = D.Propu.T_sls * thrust_SLSconv(TAS, h, D.Propu.engine.BPR, D.Propu.engine.G);

		% Equivalent plane weight [N].
		W = D.Plane.MTOW * C.g * n;

		% Drag [N].  TODO: modify when we have more precise values.
		%
		% Drag on the wings. We take only the wing induced drag as a first
		% estimation.
		CD_wing = D.Wing.CL_cr^2 / (D.Wing.e * pi * D.Wing.AR);
		D_wing  = 0.5 * rho * TAS^2 * D.Wing.S * CD_wing;
		% Drag on the fuselage. We take the parasitic drag at zero lift of the
		% entire plane, as a first estimation.
		CD_body = D.Plane.CD_0;
		D_body  = 0.5 * rho * TAS^2 * D.Wing.S * CD_body;  % TODO: check if wing surf or plane wetted surf.

		% Pitching moment [N*m].
		% This formula does not apply to our geometry.
		% C_M_0 = -0.1;
		% K_n   = 0.09;
		% C_M = C_M_0 - K_n * D.Wing.CL_max - 0.0015 * psi_max;
		C_M = D.Wing.airfoil.C_m;
		M = 0.5 * C_M * TAS^2 * rho * D.Wing.S * D.Wing.smc;

		I_theta = 3000; % TODO: to be computed with the CAD.
		c_l_alpha_plane = 6.38;
		alpha_L0_f = -0.0539;

		% Angle of attack, wing lift and HT load are coupled through three
		% equations. We solve them numerically, thanks to vpasolve().
		% TODO: change once we have the z_pos of the required components
		% (wait for the CAD).
		syms aoa L P;
		eqI   = c_l_alpha_plane * (aoa - alpha_L0_f) - L / (0.5 * C.rho_cr * D.Wing.S * (EAS^2)); %Lift of the wing
		eqII  = abs(c_g - pos(1)) * L - abs(z_cg - z_pos(8)) * T + z_cg * D_body - z_pos(1) * D_wing - abs(pos(4) - c_g) * P + M(idx) - I_theta * theta_dd_add; % Moments about COG
		eqIII = L + P - W + T(idx) * sin(aoa-alpha_L0_f); % Vertical equilibrium
		% TODO: changes names to not mess up with syms.
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