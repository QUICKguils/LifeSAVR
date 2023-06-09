% TODO:
% - The spars web should not have the same tichness as the skin.
% - We neglected the air induction system.
% - data structures in wing_geometry are a bit messy.

function struct_stress(opts)
% STRUCT_STRESS  Structural stresses and resulting design choices.
%
% This function computes the structural stresses evolving on the
% selected wing and fuselage cross sections, for all the critical points
% of the flight envelope. The resulting design choices are also yielded.
%
% Arguments:
%	opts: char {'p', 'w'}, optional
%	  'p' -> Enable plots creation.
%	  'w' -> Write plotting data in external file.
%     Default is 'p'.
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

% Utilities.
addpath(genpath(fullfile(file_dir, "../Utils")));

%% Options setting

% Option default: generate the plot.
if ~nargin
	opts = 'p';
end

%% Main solve

% 1. General data.

% Material properties.
mu        = 29.4e9;                % Shear modulus [Pa].
sigma_y   = 510e6;                 % 0.1 % proof stress [Pa].
tau_y     = sigma_y/sqrt(3);       % Shear strength [Pa] 
sigma_max = sigma_y / C.K_safety;  % Max. axial stress allowed [Pa].
tau_max   = tau_y   / C.K_safety;  % Max. shear stress allowed [Pa].

% 2. Fuselage computations.

% Fuselage data.
ns.fus = 24;  % Number of stringers.

% This table contains the relevant stresses and resulting dimensions in
% the selected fuselage cross sections, for al the CP.
FusStresses = table(...
	zeros(5, 1), zeros(5, 1), zeros(5, 1), zeros(5, ns.fus), zeros(5, 1), zeros(5, ns.fus), zeros(5, 1), ...
	'VariableNames', {'x', 'n', 'EAS', 'B_sigma_xx', 'B_min', 'q', 't_min'});

% Get the fuselage geometry parameters.
FusGeo = fus_geometry;  % TODO: should technically loop for all cs.

% Compute stresses for all the cs and CP.
for row = 1:height(D.FusLoads)
	fl = D.FusLoads(row, :);
	FusStresses(row, :) = fuselage_stresses(FusGeo, fl);
end

% Relevant fuselage design data.
FusDesign.MinBoomArea      = max(FusStresses.B_min);
FusDesign.MinSkinThickness = max(FusStresses.t_min);

% Save relevant data structures in data.mat.
save(fullfile(file_dir, "../data.mat"), "FusStresses", 'FusDesign', "-append");

% 3. Wing computations.

% Wing data.
ns.P2 = 20;  % Number of stringers in panel 2.
ns.P3 = 20;  % Number of stringers in panel 3.
ns.P4 = 20;  % Number of stringers in panel 4.
ns.wing = ns.P2 + ns.P3 + ns.P4;
t_ratio = 2.5;  % Thickness ratio between spars web and skins.  % TODO: not accurate method.

% This table contains the relevant stresses and resulting dimensions in
% the selected wing cross sections, for al the CP.
WingStresses = table(...
	zeros(5, 1), zeros(5, 1), zeros(5, 1), zeros(5, ns.wing), zeros(5, 1), zeros(5, ns.wing), zeros(5, 1), ...
	'VariableNames', {'y', 'n', 'EAS', 'B_sigma_yy', 'B_min', 'q', 't_min'});

% Get the wing geometry parameters.
WingGeo = wing_geometry;  % TODO: should technically loop for all cs.

% Compute stresses for all the cs and CP.
for row = 1:height(D.WingLoads)
	wl = D.WingLoads(row, :);
	WingStresses(row, :) = wing_stresses(WingGeo, wl);
end

% Relevant wing design data.
WingDesign.MinBoomArea      = max(WingStresses.B_min);
WingDesign.MinSkinThickness = max(WingStresses.t_min);

% Save relevant data structures in data.mat.
save(fullfile(file_dir, "../data.mat"), "WingGeo", "WingStresses", 'WingDesign', "-append");

% 4. Plot and save data.

if contains(opts, 'p')
	plot_write(FusGeo, WingGeo);
end

%% Fuselage

	function stresses = fuselage_stresses(fg, fl)
		% FUSELAGE_STRESSES  Stresses in the fuselage.
		%
		% This function implements the methodology and formulae that can be
		% found at: Structures>lessson 6>slides 41 to 63.

		% 1. Axial stress sigma_xx [Pa] and boom area B [m²].
		%
		% The axial stresses in the booms are induced by the bending moments
		% fl.My and fl.Mz. These moments generate an equivalent normal load
		% B*sigma_xx in each boom, called B_sigma_xx.
		%
		% B_sigma_xx is calculated for each booms. B is then chosen so that the
		% maximum sigma_xx does not exceed sigma_max.

		% Normal load induced by bending moments [N].
		% Obtained thanks to the Navier bending equation.
		B_sigma_xx = fl.My.*fg.z./fg.I.yy2B - fl.Mz.*fg.y./fg.I.zz2B;

		% Minimum area of the stringers [m²].
		B_min = max(abs(B_sigma_xx)) / sigma_max;

		% 2. Shear stress tau [Pa] and skin thickness t [m].
		%
		% The shear stresses in the booms are induced by the shear loads fl.Ty
		% and fl.Tz, and the torsion moment fl.Mx. These loads generate a shear
		% flow tau*t in each boom, called q.
		%
		% q is calculated along the fuselage perimeter. t is then choseDn so that
		% the maximum tau does not exceed tau_max.

		% Shear flow induced by shear loads [N/m].
		%
		% Preallocation.
		q_Ty = zeros(size(fg.y));
		q_Tz = zeros(size(fg.z));
		% First value.  TODO: these formulae are not correct.
		q_Ty(1) = -0.5 * fl.Ty/fg.I.zz2B * fg.y(1);
		q_Tz(1) = -0.5 * fl.Tz/fg.I.yy2B * fg.z(1);
		% Iteratively compute q along the fuselage perimeter.
		for i = 2:ns.fus
			q_Ty(i) = q_Ty(i-1) - fl.Ty/fg.I.zz2B * fg.y(i);
			q_Tz(i) = q_Tz(i-1) - fl.Tz/fg.I.yy2B * fg.z(i);
		end

		% Shear flow induced by torsion moment [N/m].
		q_Mx = fl.Mx / (2 * fg.A);

		% Total shear flow [N/m].
		q_tot = q_Ty + q_Tz + q_Mx;

		% Minimum thickness of the skin [m].
		t_min = max(abs(q_tot)) / tau_max;

		% Return the computed stresses and designs.
		stresses = table(fl.x, fl.n, fl.EAS, B_sigma_xx, B_min, q_tot, t_min);
	end

	function FusGeo = fus_geometry
		% FUS_GEOMETRY  Build fuselage geometry elements.

		% We choose an constant angular distribution of the booms.
		% It gives a slightly denser distribution at the top and at the
		% bottom. This is rather good: the moment of area is the weakest
		% in the vertical axis, while the bending moment is the largest
		% -> good to reinforce these locations.
		theta = 0 : 2*pi/ns.fus : 2*pi*(ns.fus-1)/ns.fus;

		% Fuselage elliptic cross section.
		a  = D.Fus.a;  % Fuselage semi major axis, in Y-direction [m].
		b  = D.Fus.b;  % Fuselage semi minor axis, in Z-direction [m].

		% Distances between the stringers and the ellipse center [m].
		r = @(theta) a*b ./ sqrt((b.*cos(theta)).^2 + (a.*sin(theta)).^2);
		% Radial distances of the stringers [m].
		r_stringers = r(theta);
		% YZ coordinates of the stringers [m].
		y = r_stringers.*cos(theta);
		z = r_stringers.*sin(theta);

		% Moment of area per unit boom area [m²].
		% Note that product moment of area Iyz is null, by symmetry.
		I.yy2B = sum(z.^2);
		I.zz2B = sum(y.^2);

		% Build return data structure.
		FusGeo.a    = a;
		FusGeo.b    = b;
		FusGeo.r    = r;
		FusGeo.y    = y;
		FusGeo.z    = z;
		FusGeo.I    = I;
		FusGeo.A    = pi*a*b;
		FusGeo.CG.x = 0;
		FusGeo.CG.z = 0;
	end

%% Wing

	function stresses = wing_stresses(wg, wl)
		% WING_STRESSES  Stresses in the wings.
		%
		% This function implements the methodology and formulae that can be
		% found at: Aircraft Structures>lessson 5>slides 15 to 52.

		% 0. Design data and preliminary computations.

		% Set frame of reference at the centroid.  % TODO: not super clean
		x_af  = wg.S.AF.x - wg.CG.x;
		z_af  = wg.S.AF.z - wg.CG.z;
		x_P2  = wg.S.P2.x - wg.CG.x;
		z_P2  = wg.S.P2.z - wg.CG.z;
		x_P3  = wg.S.P3.x - wg.CG.x;
		z_P3  = wg.S.P3.z - wg.CG.z;
		x_P4  = wg.S.P4.x - wg.CG.x;
		z_P4  = wg.S.P4.z - wg.CG.z;
		Ixx2B = wg.I.xx2B - wg.CG.z^2;
		Izz2B = wg.I.zz2B - wg.CG.x^2;
		Ixz2B = wg.I.xz2B - wg.CG.x*wg.CG.z;

		% 1. Axial stress sigma_yy [Pa] and boom area B [m²].
		%
		% The axial stresses in the booms are induced by the bending moments
		% wl.Mx and wl.Mz. These moments generate an equivalent normal load
		% B*sigma_xx in each boom, called B_sigma_xx.
		%
		% B_sigma_yy is calculated for each booms. B is then chosen so that the
		% maximum sigma_xx does not exceed sigma_max.

		% Normal load induced by bending moments [N].
		% Structures>lesson 5>slide 17.
		B_sigma_yy = ...
			(  (Izz2B*wl.Mx + Ixz2B*wl.Mz) .* z_af   ...
			  -(Ixz2B*wl.Mx + Ixx2B*wl.Mz) .* x_af ) ...
			/(Ixx2B*Izz2B - Ixz2B^2);

		% Minimum area of the stringers [m²].
		B_min = max(abs(B_sigma_yy)) / sigma_max;

		% 2. Shear stress tau [Pa] and skin thickness t [m].
		%
		% The shear stresses in the booms are induced by the shear loads wl.Tx
		% and wl.Tz, and the torsion moment wl.My. These loads generate a shear
		% flow tau*t in each boom, called q.
		%
		% q is calculated along the wing perimeter. t is then chosen so that
		% the maximum tau does not exceed tau_max.

		% 2.1. Shear flow induced by torsion moment [N/m].
		%
		% We compute them by using the compatibility of twist rate.
		% Structures>lesson 5>slides 18 to 25.

		% The torsion shear flows in the two cells are coupled with the
		% twist rate of the profile through three equations:
		% - twist rate equilibrium in each cell,
		% - torsion moment equilibrium of the profile.

		% Define symbolic variables.
		tr = sym('tr');  % Twist rate [rad/s].
		q1 = sym('q1');  % Shear flow in cell 1 [N/m].
		q2 = sym('q2');  % Shear flow in cell 2 [N/m].
		q3 = sym('q3');  % Shear flow in cell 3 [N/m].

		% Define system of equations.
		eqns = [ ...
			tr == 1/(2*wg.A(1)*mu) * (wg.L(3)*q1 + wg.L(6)*(q1-q2)), ...
			tr == 1/(2*wg.A(2)*mu) * ((wg.L(2)+wg.L(4))*q2 + wg.L(6)*(q2-q1) + wg.L(7)*(q2-q3)), ...
			tr == 1/(2*wg.A(3)*mu) * ((wg.L(1)+wg.L(5))*q3 + wg.L(7)*(q3-q2)), ...
			wl.My == - 2 * (wg.A(1)*q1 + wg.A(2)*q2 + wg.A(3)*q3);
		];

		% Solve the system.
		res = vpasolve(eqns, [tr, q1, q2, q3]);

		% Unpack the results of interest.
		q_My(1) = double(res.q1);
		q_My(2) = double(res.q2);
		q_My(3) = double(res.q3);

		% 2.2. Shear flow induced by shear loads [N/m].
		%
		% Structures>lesson 5>slides 26 to 52.

		% 2.2.1. Taper and sweep effect.
		% Structures>lesson 5>slide 40.

		% x and z derivatives of the stringers along y.
		[dxdy, dzdy] = derive_stringers(x_af, z_af);

		% Web shear loads.
		% This gives the "corrected" shear loads exerting on a web
		% subtended by the stringers of indexes i1_af and i2_af, in the
		% index datum of the airfoil.
		Tx_w = @(i2_af, i1_af) wl.Tx - B_sigma_yy(i2_af)*dxdy(i2_af) - B_sigma_yy(i1_af)*dxdy(i1_af);
		Tz_w = @(i2_af, i1_af) wl.Tz - B_sigma_yy(i2_af)*dzdy(i2_af) - B_sigma_yy(i1_af)*dzdy(i1_af);

		% 2.2.2. Open shear flows.
		%
		% We cut the cells at the top left of each spar booms.

		% Incremental shear flow in skins.
		qo = @(i2_af, i1_af) ...
			- (Izz2B*Tz_w(i2_af, i1_af)-Ixz2B*Tx_w(i2_af, i1_af))/(Ixx2B*Izz2B-Ixz2B^2) * z_af(i1_af) ...
			- (Ixx2B*Tx_w(i2_af, i1_af)-Ixz2B*Tz_w(i2_af, i1_af))/(Ixx2B*Izz2B-Ixz2B^2) * x_af(i1_af);

		% Open shear flow in panel 3.
		qo_P3 = zeros(1, ns.P3-1);
		qo_P3(1) = 0;  % Because of the cut at this place.
		for i = 2:ns.P3-1
			qo_P3(i) = qo_P3(i-1) + qo(ns.P2+(i+1), ns.P2+i);
		end

		% Open shear flow in panel 2.
		qo_P2 = zeros(1, ns.P2-1);
		qo_P2(1) = 0;  % Because of the cut at this place.
		for i = 2:ns.P2-1
			qo_P2(i) = qo_P2(i-1) + qo(i+1, i);
		end

		% Open shear flow in panel 6.
		qo_P6 = qo_P2(end) + 2 * qo(ns.P2+ns.P3, ns.P2);  % Continues from panel 2.

		% Open shear flow in panel 7.
		qo_P7 = qo(ns.wing, 1);  % Just the first stringer, next to the cut.
	
		% Open shear flow in panel 4.
		qo_P4 = zeros(1, ns.P4-1);
		qo_P4(1) = qo_P7 + qo(ns.wing-1, ns.wing);  % Continues from panel 7.
		for i = 2:ns.P4-1
			qo_P4(i) = qo_P4(i-1) + qo(ns.wing-i, ns.wing-(i-1));
		end

		% 2.2.3. Closed shear flows at cut (correction).
		%
		% We compute them by using the compatibility of twist rate.
		%
		% The closed shear flows in the two cells are coupled with the
		% twist rate of the profile through three equations:
		% - twist rate equilibrium in each cell,
		% - torsion moment equilibrium of the profile.

		% Define symbolic variables.
		tr = sym('tr');  % Twist rate [rad/s].
		q1 = sym('q1');  % Closed shear flow in cell 1 [N/m].
		q2 = sym('q2');  % Closed shear flow in cell 2 [N/m].
		q3 = sym('q3');  % Closed shear flow in cell 3 [N/m].

		% Average distance between two succesive stringers.  % PERF: slow.
		d = @(x, z) (vecnorm([x(2:end); z(2:end)]) - vecnorm([x(1:end-1); z(1:end-1)]))/2;

		% Define system of equations.
		eqns = [ ...
			tr == 1/(2*wg.A(1)*mu) * ( ...
				wg.L(3)*q1 + wg.L(6)*(q1-q2) ...
				+ sum(qo_P3) * wg.L(3)/(ns.P3-1) - qo_P6 * wg.L(6)), ...
			tr == 1/(2*wg.A(2)*mu) * ( ...
				(wg.L(2)+wg.L(4))*q2 + wg.L(6)*(q2-q1) + wg.L(7)*q2 ...
				+ sum(qo_P2) * wg.L(2)/(ns.P2-1) - sum(qo_P4) * wg.L(4)/(ns.P4-1) + qo_P6 * wg.L(6) - qo_P7 * wg.L(7)), ...
			tr == 1/(2*wg.A(3)*mu) * ( ...
				(wg.L(1)+wg.L(5))*q3 + wg.L(7)*(q3-q2) + qo_P7 * wg.L(7)), ...
			wl.My == ...
				  sum(d(x_P2, z_P2).*qo_P2) * wg.L(2)/(ns.P2-1) ...
				+ sum(d(x_P3, z_P3).*qo_P3) * wg.L(3)/(ns.P3-1) ...
				- sum(d(x_P4, z_P4).*qo_P4) * wg.L(4)/(ns.P4-1) ...
				+ 2 * qo_P6 * wg.Ah(6) ...
				- 2 * qo_P7 * wg.Ah(7) ...
				+ 2 * (wg.Ah(3) - wg.Ah(6)) * q1 ...
				+ 2 * (wg.Ah(2) + wg.Ah(6) - wg.Ah(4) - wg.Ah(7)) * q2 ...
				+ 2 * (wg.Ah(1) + wg.Ah(7) + wg.Ah(5)) * q3 ...
				+ sum(x_af .* B_sigma_yy .* dzdy) ...
				- sum(z_af .* B_sigma_yy .* dxdy) ...
		];

		% Solve the system.
		res = vpasolve(eqns, [tr, q1, q2, q3]);

		% Unpack the results of interest.
		qc(1) = double(res.q1);
		qc(2) = double(res.q2);
		qc(3) = double(res.q3);

		% 2.3. Total shear flow and miminum skin thickness.

		% Total shear flows [N/m].
		q_P2  =   q_My(2)         + qo_P2 + qc(2);
		q_P3  =   q_My(1)         + qo_P3 + qc(1);
		q_P4  = - q_My(2)         + qo_P4 - qc(2);
		q_P6  = - q_My(1)+q_My(2) + qo_P6 - qc(1)+qc(2);
		q_P7  = - q_My(2)+q_My(3) + qo_P7 - qc(2)+qc(3);
		q_P15 =   q_My(3)                 + qc(3);

		% Gather them in one vector.
		q_tot = [q_P2, q_P3, q_P4, q_P6, q_P7, q_P15];
		% Discard the spars web, as they are not technically "skins".
		q_tot_skin = [q_P2, q_P3, q_P4, q_P15] / t_ratio;

		% Minimum thickness of the skin [m].
		t_min = max(abs(q_tot_skin)) / tau_max;

		% Return the computed stresses and designs.
		stresses = table(wl.y, wl.n, wl.EAS, B_sigma_yy, B_min, q_tot, t_min);
	end

	function WingGeo = wing_geometry
		% WING_GEOMETRY  Build wing geometry elements.

		% Load the file containing the X and Z coordinates of the airfoil
		% profile, starting from the leading edge and cycling clockwise.
		airfoil = load(fullfile(file_dir, "sc2_0412.dat"));
		% Extract XZ coordinates of the perimeter [m].
		x_af = D.Wing.c_root * airfoil(:,1);
		z_af = D.Wing.c_root * airfoil(:,2);
		% Parametric curve of the airfoil [m].
		AF = new_pc(x_af, z_af);

		% Spars X coord.
		x_spar1 = D.Wing.c_root * (1 - D.Wing.flaps.chord_ratio);
		x_spar2 = D.Wing.c_root / 4;  % TODO replace with the aerodynamic center.

		% Nodes of the profile.
		% We call "nodes" the intersections of the spars with the airfoil.
		N{1}.t = [fzero(@(t) AF.x(t) - x_spar1, [0, 0.5])];
		N{2}.t = [fzero(@(t) AF.x(t) - x_spar2, [0, 0.5])];
		N{3}.t = [fzero(@(t) AF.x(t) - x_spar2, [0.5, 1])];
		N{4}.t = [fzero(@(t) AF.x(t) - x_spar1, [0.5, 1])];
		N{1}.x = x_spar1;
		N{2}.x = x_spar2;
		N{3}.x = x_spar2;
		N{4}.x = x_spar1;
		N{1}.z = AF.z(N{1}.t);
		N{2}.z = AF.z(N{2}.t);
		N{3}.z = AF.z(N{3}.t);
		N{4}.z = AF.z(N{4}.t);

		% Panels of the profile.
		% We call "panels" the walls perpendicular to the airfoil cross section.
		P{1} = extract_pc(AF, 0,      N{1}.t);
		P{2} = extract_pc(AF, N{1}.t, N{2}.t);
		P{3} = extract_pc(AF, N{2}.t, N{3}.t);
		P{4} = extract_pc(AF, N{3}.t, N{4}.t);
		P{5} = extract_pc(AF, N{4}.t, 1);
		P{6} = new_pc([N{2}.x, N{3}.x], [N{2}.z, N{3}.z]);
		P{7} = new_pc([N{1}.x, N{4}.x], [N{1}.z, N{4}.z]);
		% Panel lengths.
		L = cellfun(@pc_length, P);

		% Cells of the profile.
		% We call "cells" the two closed sections of the airfoil profile.
		Cell{1} = join_pc(P{3}, reverse_pc(P{6}));
		Cell{2} = join_pc(join_pc(P{2}, P{6}), join_pc(P{4}, reverse_pc(P{7})));
		Cell{3} = join_pc(join_pc(P{2}, P{3}), join_pc(P{4}, reverse_pc(P{7})));

		% Stringers of the profile.
		% We call "stringers"... seriously, take a look at the course bruh.
		%
		% On panel 2.
		S.P2.t = linspace(0, 1, ns.P2);
		S.P2.x = P{2}.x(S.P2.t);
		S.P2.z = P{2}.z(S.P2.t);
		% On panel 3.
		S.P3.t = linspace(0, 1, ns.P3);
		S.P3.x = P{3}.x(S.P3.t);
		S.P3.z = P{3}.z(S.P3.t);
		% On panel 4.
		S.P4.t = linspace(0, 1, ns.P4);
		S.P4.x = P{4}.x(S.P4.t);
		S.P4.z = P{4}.z(S.P4.t);
		% On profile.
		S.AF.t = [S.P2.t, S.P3.t, S.P4.t];
		S.AF.x = [S.P2.x, S.P3.x, S.P4.x];
		S.AF.z = [S.P2.z, S.P3.z, S.P4.z];

		% Centroid of the profile [m].
		x_CG = sum(S.AF.x) / ns.wing;
		z_CG = sum(S.AF.z) / ns.wing;

		% Moment of area per unit boom area of the profile [m²].
		I.xx2B = sum(S.AF.z.^2);
		I.zz2B = sum(S.AF.x.^2);
		I.xz2B = sum(S.AF.x .* S.AF.z);

		% Swept area of the panels, from the CG.
		Ah = cellfun(@(pc) swept_area(translate_pc(pc, x_CG, z_CG)), P);
		% Cell areas.
		A(1) = Ah(3) - Ah(6);
		A(2) = Ah(2) + Ah(6) + Ah(4) + Ah(7);
		A(3) = Ah(1) + Ah(5) - Ah(7);
		% NOTE: computing cell areas with `cellfun(@enclosed_area, Cell)` would
		% be more elegant, but it leads to perf issues, and I don't know why.

		% Build return data structure.
		WingGeo.N    = N;
		WingGeo.P    = P;
		WingGeo.Ah   = Ah;
		WingGeo.L    = L;
		WingGeo.C    = Cell;
		WingGeo.A    = A;
		WingGeo.S    = S;
		WingGeo.I    = I;
		WingGeo.CG.x = x_CG;
		WingGeo.CG.z = z_CG;
		WingGeo.x_af = x_af;
		WingGeo.z_af = z_af;
	end

	function pc = new_pc(Xs, Zs)
		% NEW_PC  Create new parametric curve.
		%
		% Arguments:
		%   Xs (1xn double) -- X coord. sample of the curve.
		%   Ys (1xn double) -- Z coord. sample of the curve.
		% Return:
		%   pc (function_handle) -- Parametric curve.
		%     Argument:
		%       t -- Parameter, ranging from 0 to 1.
		%     Return:
		%       s -- Position vector.
		%       x -- X component of s.
		%       z -- Z component of s.
		pc.x = @(t) interp1(linspace(0, 1, numel(Xs)), Xs, t, 'linear', 'extrap');
		pc.z = @(t) interp1(linspace(0, 1, numel(Zs)), Zs, t, 'linear', 'extrap');
		pc.s = @(t) [pc.x(t), pc.z(t)];
	end

	function pc_reversed = reverse_pc(pc)
		% INVERT_PC  Reverse the path direction of the parametric curve.
		pc_reversed.s = @(t) pc.s(1-t);
		pc_reversed.x = @(t) pc.x(1-t);
		pc_reversed.z = @(t) pc.z(1-t);
	end

	function pc_translated = translate_pc(pc, x, z)
		% TRANSLATE_PC  Translate the parametric curve.
		pc_translated.x = @(t) pc.x(t) - x;
		pc_translated.z = @(t) pc.z(t) - z;
		pc_translated.s = @(t) [pc_translated.x(t), pc_translated.z(t)];
	end

	function pc_extracted = extract_pc(base_pc, ti, tf)
		% EXTRACT_PC  Extract a parametric curve from an existing one.
		pc_extracted.x = @(t) base_pc.x(ti + t*(tf-ti));
		pc_extracted.z = @(t) base_pc.z(ti + t*(tf-ti));
		pc_extracted.s = @(t) base_pc.s(ti + t*(tf-ti));
	end

	function pc_joined = join_pc(pc1, pc2)
		% JOIN_PC  Join two parametric curves.
		pc_joined.x = @(t) pc1.x(2*t) * (t<=0.5) + pc2.x(2*(t-0.5)) * (t>0.5);
		pc_joined.z = @(t) pc1.z(2*t) * (t<=0.5) + pc2.z(2*(t-0.5)) * (t>0.5);
		pc_joined.s = @(t) pc1.s(2*t) * (t<=0.5) + pc2.s(2*(t-0.5)) * (t>0.5);
	end

	function pc_derived = derive_pc(pc, dt)
		% DERIVE_PC  Approximate derivative of a parametric curve.

		% Default value for dt.
		if nargin == 1; dt = 1e-5; end

		pc_derived.x = @(t) (pc.x(t+dt) - pc.x(t)) ./ dt;
		pc_derived.z = @(t) (pc.z(t+dt) - pc.z(t)) ./ dt;
		pc_derived.s = @(t) (pc.s(t+dt) - pc.s(t)) ./ dt;
	end

	function length = pc_length(pc, ti, tf)
		% PC_LENGTH  Arc length of a parametric curve.

		% Default: calculate the length of the entire curve.
		if nargin == 1; ti = 0; tf = 1; end

		pc_derived = derive_pc(pc);
		length = integral(@(t) norm(pc_derived.s(t)), ti, tf, 'ArrayValued', true);
	end

	function area = swept_area(pc, ti, tf)  % PERF: too slow
		% SWEPT_AREA  Arc area enclosed by a parametric curve.
		%
		% Green's formula: A = 0.5 * int_ti^tf((x(t)y'(t)-x'(t)y(t))dt.

		% Default: calculate the area enclosed by the entire curve.
		if nargin == 1; ti = 0; tf = 1; end

		dCdt = derive_pc(pc);
		integrand = @(t) pc.x(t) * dCdt.z(t) - dCdt.x(t) * pc.z(t);
		area = 0.5 * integral(integrand, ti, tf, 'ArrayValued', true);
		area = abs(area);
	end

	function [dxdy, dzdy] = derive_stringers(x_cs, z_cs)
		% DERIVE_STRINGERS  Spanwise x and z derivatives of the stringers.

		% Length of the spars [m].
		L = (D.Wing.span/2 - wl.y) / cosd(D.Wing.sweep);

		% X-translation of the profile, due to sweep [m].
		% tr = (D.Wing.span/2 - wl.y) * tand(D.Wing.sweep);

		% Scaling of the profile, due to taper.
		sc = D.Wing.taper;

		% Apply affine transformation: scaling and translation.
		% x_tip = x_cs .* sc + tr;  % TODO: not sure if we should take sweep into account.
		x_tip = x_cs .* sc;
		z_tip = z_cs .* sc;

		% Debug print: uncomment to see tip airfoil profile.
		% hold on; plot(x_tip, z_tip, 'Marker','+'); hold off;

		% Get the spanwise derivatives.
		dxdy = (x_tip - x_cs) ./ L;
		dzdy = (z_tip - z_cs) ./ L;
	end

	function plot_write(fg, wg)
		% PLOT_WRITE  Plot the wing's geometry elements of the profile.

		% Instantiate a figure object.
		figure('WindowStyle','docked');
		hold on;

		% Panels.
		cellfun(@(p) fplot(p.x, p.z, [0, 1]), wg.P);
		% Nodes.
		cellfun(@(n) plot(n.x, n.z, 'o'), wg.N);
		% Stringers.
		plot(wg.S.P2.x, wg.S.P2.z, '*');
		plot(wg.S.P3.x, wg.S.P3.z, '*');
		plot(wg.S.P4.x, wg.S.P4.z, '*');
		% COG.
		plot(wg.CG.x, wg.CG.z, 'x');

		% Dress the plot.
		title('Geometry of the wing section');
		grid;
		axis equal;
		hold off;


		% Write plotting data in external file, if desired.
		if contains(opts, 'w')
			write_ext();
		end

		function write_ext()
			% WRITE_EXT	 Write plotting data in external file.

			% Generate the filename of the object to save.
			gen_fname = @(obj) fullfile(file_dir, strcat("../Plot/", obj, ".dat"));

			% Fuselage contour.
			writematrix([linspace(0, 2*pi, 100)', arrayfun(fg.r, linspace(0, 2*pi, 100))'], gen_fname("FusContour"));
			% Fuselage stringers.
			writematrix([fg.y', fg.z'], gen_fname("FusStringers"));
			% Wing contour.
			writematrix([wg.x_af, wg.z_af], gen_fname("WingContour"));
			% Wing stringers.
			writematrix([wg.S.AF.x', wg.S.AF.z'], gen_fname("WingStringers"));
		end
	end
end