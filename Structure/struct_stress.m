% TODO:
% - Some values are guessed. Wait for more precise values.
% - Lots of hardcoded values.
% - Lots of computation are redudant and should not be done for each
%   Loads table lines. For examle, the stringers coordinates.
% - We neglected the air induction system.
% - the pc structure is not very clean. Fins a way to properly extract x
%   and z from s, without explicitely define their expressions.

function struct_stress
% STRUCT_STRESS  Structural stresses and resulting design choices.
%
% This function computes the structural stresses evolving on the
% selected wing and fuselage cross sections, for all the critical points
% of the flight envelope. The resulting design choices are also yielded.
%
% This function implements the methodology and formulae that can be
% found at: Aircraft Structures>lessson 5>slides 41 to 63.
%
% Save:
%   FusStresses  (table)  -- Stresses in the fuselage cross sections.
%   WingStresses (table)  -- Stresses in the wings    cross sections.
%   FusDesign    (struct) -- Resulting design choices for the fuselage.
%   WingDesign   (struct) -- Resulting design choices for the wing.

%% Imports

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Constants and aircraft data.
C = load(fullfile(file_dir, "../constants.mat"));
D = load(fullfile(file_dir, "../data.mat"));

% Aliases.
a  = D.Fus.a;      % Fuselage semi major axis, in Y-direction [m].
b  = D.Fus.b;      % Fuselage semi minor axis, in Z-direction [m].
FL = D.FusLoads;   % MNT in the selected fuselage cs, for each CP of the FE.
% WL = D.WingLoads;  % MNT in the selected wing     cs, for each CP of the FE.

% Utilities.
addpath(genpath(fullfile(file_dir, "../Utils")));

%% Main solve

% 1. General data.

% Material properties.  TODO: update when matsel is done.
sigma_y   = 280e6;                 % 0.1 % proof stress [Pa].
tau_y     = sigma_y/sqrt(3);       % Shear strength [Pa]
sigma_max = sigma_y / C.K_safety;  % Max. axial stress allowed [Pa].
tau_max   = tau_y   / C.K_safety;  % Max. shear stress allowed [Pa].

% 2. Fuselage computations.

% Fuselage data.
ns_fus = 24;  % Number of stringers.  TODO: see later if it yield decent booms area.

% This table contains the relevant stresses and resulting dimensions in
% the selected fuselage cross sections, for al the CP.
FusStresses = table(...
	zeros(5, 1), zeros(5, 1), zeros(5, 1), zeros(5, 24), zeros(5, 1), zeros(5, 24), zeros(5, 1), ...
	'VariableNames', {'x', 'n', 'EAS', 'B_sigma_xx', 'B_min', 'q', 't_min'});

% Compute stresses for all the cs and CP.
for row = 1:height(FL)
	fl = FL(row, :);
	FusStresses(row, :) = fuselage_stresses(fl);
end

FusDesign.MinBoomArea     = max(FusStresses.B_min);
FusDesign.MinSkinTickness = max(FusStresses.t_min);

% Save relevant data structures in data.mat.
save(fullfile(file_dir, "../data.mat"), "FusStresses", 'FusDesign', "-append");

% 3. Wing computations.

% % Wing data.
% ns_wing = 39;  % Number of stringers.  TODO: see later if it yield decent booms area.
%
% % This table contains the relevant stresses and resulting dimensions in
% % the selected wing cross sections, for al the CP.
% WingStresses = table(...
% 	zeros(5, 1), zeros(5, 1), zeros(5, 1), zeros(5, 24), zeros(5, 1), zeros(5, 24), zeros(5, 1), ...
% 	'VariableNames', {'x', 'n', 'EAS', 'B_sigma_xx', 'B_min', 'q', 't_min'});
%
% % Compute stresses for all the cs and CP.
% for row = 1:height(WL)
% 	wl = WL(row, :);
% 	WingStresses(row, :) = wing_stresses(wl);
% end
%
% WingDesign = {};
wings_stresses();

% 4. Save relevant data structures in data.mat.

save(fullfile(file_dir, "../data.mat"), "FusStresses",  'FusDesign',  "-append");
% save(fullfile(file_dir, "../data.mat"), "WingStresses", 'WingDesign', "-append");

%% Fuselage

	function stresses = fuselage_stresses(fl)
		% FUSELAGE_STRESSES  Stresses in the fuselage.

		% 0. Design data and preliminary computations.

		% We choose an constant angular distribution of the booms.
		% Ca donne une distrib. légèrement plus dense au dessus et en
		% dessous. C'est pas plus mal, c'est dans ce sens que l'inertie
		% est la plus faible et que les moments de flexion sont lles plus
		% important -> bien de renforcer ces endroits.
		theta = 0 : 2*pi/ns_fus : 2*pi*(ns_fus-1)/ns_fus;

		% Distances between the stringers and the ellipse center [m].
		r = a*b ./ sqrt((b.*cos(theta)).^2 + (a.*sin(theta)).^2);
		% YZ coordinates of the stringers [m].
		y = r.*cos(theta);
		z = r.*sin(theta);
		% Uncomment line below to check the stringers distribution.
		% figure; plot(y, z, '*'); grid;

		% 1. Axial stress sigma_xx [Pa] and boom area B [m²].
		%
		% The axial stresses in the booms are induced by the bending moments
		% al.My and al.Mz. These moments generate an equivalent normal load
		% B*sigma_xx in each boom, called B_sigma_xx.
		%
		% B_sigma_xx is calculated for each booms. B is then chosen so that the
		% maximum sigma_xx does not exceed sigma_max.

		% Moment of area per unit boom area [m²].
		% Note that product moment of area Iyz is null, by symmetry.
		Iyy2B = sum(z.^2);
		Izz2B = sum(y.^2);

		% Normal load induced by bending moments [N].
		% Obtained thanks to the Navier bending equation.
		B_sigma_xx = fl.My.*z./Iyy2B - fl.Mz.*y./Izz2B;

		% Minimum area of the stringers [m²].
		B_min = max(abs(B_sigma_xx)) / sigma_max;

		% 2. Shear stress tau [Pa] and skin thickness t [m].
		%
		% The shear stresses in the booms are induced by the shear loads fl.Ty
		% and fl.Tz, and the torsion moment fl.Mx. These loads generate a shear
		% flow tau*t in each boom, called q.
		%
		% q is calculated along the fuselage perimeter. t is then chosen so that
		% the maximum tau does not exceed tau_max.

		% Shear flow induced by shear loads [N/m].
		%
		% Preallocation.
		q_Ty = zeros(size(y));
		q_Tz = zeros(size(z));
		% First value.  TODO: these formulae are not correct.
		q_Ty(1) = -0.5 * fl.Ty/Izz2B * y(1);
		q_Tz(1) = -0.5 * fl.Tz/Iyy2B * z(1);
		% Iteratively compute q along the fuselage perimeter.
		for i = 2:ns_fus
			q_Ty(i) = q_Ty(i-1) - fl.Ty/Izz2B * y(i);
			q_Tz(i) = q_Tz(i-1) - fl.Tz/Iyy2B * z(i);
		end

		% Shear flow induced by torsion moment [N/m].
		q_Mx = fl.Mx / (2 * pi*a*b);

		% Total shear flow [N/m].
		q_tot = q_Ty + q_Tz + q_Mx;

		% Minimum tickness of the skin [m].
		t_min = max(abs(q_tot)) / tau_max;

		% Rivets sizing.
		% R = - fl.Ty/Izz2B*y_i(:,1) - fl.Tz/Iyy2B*z_i(:,1);

		% Return the computed stresses and designs.
		stresses = table(fl.x, fl.n, fl.EAS, B_sigma_xx, B_min, q_tot, t_min);
	end

%% Wing

	function wings_stresses(wl)
		% WING_STRESSES  Stresses in the wings.
		%
		% NOTE:
		%   See the report, or annex drawing that explain the labelling
		%   convention.

		% TODO: lots of computations must lie in main solve.

		% Load the file containing the X and Z coordinates of the airfoil
		% profile, starting from the leading edge and cycling clockwise.
		airfoil = load(fullfile(file_dir, "sc2_0412.dat"));
		% Extract XZ coordinates of the perimeter [m].
		x_airfoil = D.Wing.c_root * airfoil(:,1);
		z_airfoil = D.Wing.c_root * airfoil(:,2);
		% Parametric curve of the airfoil [m].
		af_pc = new_pc(x_airfoil, z_airfoil);

		% XZ coordinates of the centroïd of area of the skin [m].
		% 		x_CG = sum(x_airfoil) / numel(x_airfoil);
		% 		z_CG = sum(z_airfoil) / numel(z_airfoil);
		% Upper and lower airfoil sampling.
		% 		x_af_up   = x_airfoil(1, 1           : ceil(end/2));
		% 		x_af_down = x_airfoil(1, ceil(end/2) : end);

		% Spars X coord.
		x_spar1 = D.Wing.c_root * (1 - D.Wing.flaps.chord_ratio);
		x_spar2 = D.Wing.c_root / 4;  % TODO replace with the aerodynamic center.

		% Nodes of the profile.
		% We call "nodes" the intersections of the spars with the airfoil.
		N1.t = [fzero(@(t) af_pc.x(t) - x_spar1, [0, 0.5])];
		N1.x = x_spar1;
		N1.z = af_pc.z(N1.t);
		N2.t = [fzero(@(t) af_pc.x(t) - x_spar2, [0, 0.5])];
		N2.x = x_spar2;
		N2.z = af_pc.z(N2.t);
		N3.t = [fzero(@(t) af_pc.x(t) - x_spar2, [0.5, 1])];
		N3.x = x_spar2;
		N3.z = af_pc.z(N3.t);
		N4.t = [fzero(@(t) af_pc.x(t) - x_spar1, [0.5, 1])];
		N4.x = x_spar1;
		N4.z = af_pc.z(N4.t);

		% Parametric curves of the panels.
		% We call "panels" the
		P1 = extract_pc(af_pc, 0,    N1.t);
		P2 = extract_pc(af_pc, N1.t, N2.t);
		P3 = extract_pc(af_pc, N2.t, N3.t);
		P4 = extract_pc(af_pc, N3.t, N4.t);
		P5 = extract_pc(af_pc, N4.t, 1);
		P6 = new_pc([N1.x, N4.x], [N1.z, N4.z]);
		P7 = new_pc([N2.x, N3.x], [N2.z, N3.z]);

		% Uncomment this paragraph to check the panels and nodes.
		figure('WindowStyle','docked');
		hold on;
		fplot(P1.x, P1.z, [0, 1]);
		fplot(P2.x, P2.z, [0, 1]);
		fplot(P3.x, P3.z, [0, 1]);
		fplot(P4.x, P4.z, [0, 1]);
		fplot(P5.x, P5.z, [0, 1]);
		fplot(P6.x, P6.z, [0, 1]);
		fplot(P7.x, P7.z, [0, 1]);
		plot(N1.x, N1.z, '*'); plot(N2.x, N2.z, '*');
		plot(N3.x, N3.z, '*'); plot(N4.x, N4.z, '*');
		hold off;
		title('Panels and nodes');
		grid;
		axis equal;

		% test.
		dP4 = derive_pc(P4);
		disp(dP4.s(0.8));
		disp(pc_length(P3, 0, 1));
	end

	function pc = new_pc(Xs, Zs)
		% NEW_PC  Create new parametric curve.
		%
		% Arguments:
		%   Xs (1xn double) -- X coord. sample of the curve.
		%   Ys (1xn double) -- Z coord. sample of the curve.
		% Return:
		%   pc (function_handle) -- Parametric curve.

		pc.x = @(t) interp1(linspace(0, 1, numel(Xs)), Xs, t);
		pc.z = @(t) interp1(linspace(0, 1, numel(Zs)), Zs, t);
		pc.s = @(t) [pc.x(t), pc.z(t)];
	end

	function pc_extracted = extract_pc(base_pc, ti, tf)
		% EXTRACT_PC  Extract a parametric curve from an existing one.
		pc_extracted.x = @(t) base_pc.x(ti + t*(tf-ti));
		pc_extracted.z = @(t) base_pc.z(ti + t*(tf-ti));
		pc_extracted.s = @(t) base_pc.s(ti + t*(tf-ti));
	end

	function pc_joined = join_pc(pc1, pc2)
		% JOIN_PC  Join two parametric curves.
		pc_joined.x = @(t) pc1.x(t/2) + pc2.x(t/2+0.5);
		pc_joined.z = @(t) pc1.z(t/2) + pc2.z(t/2+0.5);
		pc_joined.s = @(t) pc1.s(t/2) + pc2.s(t/2+0.5);
	end

	function pc_derived = derive_pc(pc, dt)
		% DERIVE_PC  Approximate derivative of a parametric curve.

		% Default value for dt.
		if nargin == 1
			dt = 1e-5;
		end

		pc_derived.x = @(t) (pc.x(t+dt) - pc.x(t)) ./ dt;
		pc_derived.z = @(t) (pc.z(t+dt) - pc.z(t)) ./ dt;
		pc_derived.s = @(t) (pc.s(t+dt) - pc.s(t)) ./ dt;
	end

	function length = pc_length(pc, ti, tf)
		% PC_LENGTH  Arc length of a parametric curve.
		pc_derived = derive_pc(pc);
		length = integral(@(t) norm(pc_derived.s(t)), ti, tf, 'ArrayValued', true);
	end

	function area = cell_area(cell)
		% CELL_AREA  Area enclosed by a cell.
		%
		% Green's formula: A = 0.5 * int_0^1((x(t)y'(t)-x'(t)y(t))dt.
		area = cell;
	end
end