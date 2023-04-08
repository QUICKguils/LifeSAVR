% TODO:
% - Some values are guessed. Wait for more precise values.
% - Lots of computation are redudant and should not be done for each
%   Loads table lines. For examle, the stringers coordinates.
% - We neglected the air induction system.

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
WL = D.WingLoads;  % MNT in the selected wing     cs, for each CP of the FE.

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

% % Fuselage data.
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

		% Uncomment to check the stringers distribution.
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
		% First value.
		q_Ty(1) = -0.5 * fl.Ty/Izz2B * y(1);
		q_Tz(1) = -0.5 * fl.Tz/Iyy2B * z(1);
		% Iteratively compute q along the fuselage perimeter.
		for i = 2:ns_fus
			q_Ty(i) = q_Ty(i-1) - fl.Ty/Izz2B * y(i);
			q_Tz(i) = q_Tz(i-1) - fl.Tz/Iyy2B * z(i);
		end

		% Area of the fuselage cross section [m²].
		area = pi * a * b; 

		% Shear flow induced by torsion moment [N/m].
		q_Mx = fl.Mx / (2*area);

		% Total shear flow [N/m].
		q_tot = q_Ty + q_Tz + q_Mx;

		% Minimum tickness of the skin [m].
		t_min = max(abs(q_tot)) / tau_max;

		% Rivets sizing.
		% R_aa = (-Ty_aa/I2B_zz_aa*y_i(:,1))-(Tz_aa/I2B_yy_aa*z_i(:,1));
		% R_bb = (-Ty_bb/I2B_zz_bb*y_i(:,3))-(Tz_bb/I2B_yy_bb*z_i(:,3));
		% R_cc = (-Ty_cc/I2B_zz_cc*y_i(:,2))-(Tz_cc/I2B_yy_cc*z_i(:,2));
		% R_max = max(abs([R_aa R_bb R_cc]));

		% Return the computed stresses and designs.
		stresses = table(fl.x, fl.n, fl.EAS, B_sigma_xx, B_min, q_tot, t_min);
	end

%% Wing

	function wings_stresses(wl)
		% WING_STRESSES  Stresses in the wings.

% 		%root
% 
% 		% Load the file containing the X and Y coordinates of the airfoil
% 		% profile, starting from the leading edge and cycling clockwise.
% 		geom = load(fullfile(file_dir, "../constants.mat"));
% % 
% % 		x_l_tot = wing_l(:,1);
% % 		z_l_tot = wing_l(:,2);
% % 		x_u_tot = wing_u(:,1);
% % 		z_u_tot = wing_u(:,2);
% 
% 
% 		X_centroid = 0.4125;
% 		Z_centroid = 0.0028;
% 
% 		%spar positions
% 
% 		x_pos_last_spar = (1 - aileron_chord_ratio);
% 		x_pos_first_spar =  X_ac/c_bar_w;
% 
% 		index_last_spar_u = find(abs(x_u_tot-x_pos_last_spar)==0);
% 		index_last_spar_l = find(abs(x_l_tot-x_pos_last_spar)==0);
% 
% 		index_first_spar_u = find(abs(x_u_tot-x_pos_first_spar)<0.01);
% 		index_first_spar_l = find(abs(x_l_tot-x_pos_first_spar)<0.01);
% 
% 		x_l = x_l_tot(1:index_last_spar_l);
% 		z_l = z_l_tot(1:index_last_spar_l);
% 		x_u = x_u_tot(1:index_last_spar_u);
% 		z_u = z_u_tot(1:index_last_spar_u);
% 
% 
% 		%lower length spacing left part
% 
% 		for i = 2:index_last_spar_l
% 			l_interval_1(i) = sqrt((x_l(i)-x_l(i-1))^2 + (z_l(i)-z_l(i-1))^2);
% 		end
% 		perimeter_norm_1 = sum(l_interval_1);
% 
% 		for i = 2:index_last_spar_u
% 			l_interval_2(i) = sqrt((x_u(i)-x_u(i-1))^2 + (z_u(i)-z_u(i-1))^2);
% 		end
% 		perimeter_norm_2 = sum(l_interval_2);
% 
% 		perimeter_norm_tot = perimeter_norm_1 + perimeter_norm_2;
% 
% 
% 
% 		%% perimeters
% 		%perimeter 1
% 		for i = index_first_spar_u:index_last_spar_u
% 			l_interval_u(i) = sqrt((x_u(i)-x_u(i-1))^2 + (z_u(i)-z_u(i-1))^2);
% 		end
% 		perimeter_norm_u = sum(l_interval_u);
% 
% 		%perimeter 2
% 		for i = index_first_spar_l:index_last_spar_l
% 			l_interval_l(i) = sqrt((x_l(i)-x_l(i-1))^2 + (z_l(i)-z_l(i-1))^2);
% 		end
% 		perimeter_norm_l = sum(l_interval_l);
% 
% 		%perimeter 3
% 		for i = 2:index_first_spar_l
% 			l_interval_left_1(i) = sqrt((x_l(i)-x_l(i-1))^2 + (z_l(i)-z_l(i-1))^2);
% 		end
% 		perimeter_norm_left_1 = sum(l_interval_left_1);
% 
% 		for i = 2:index_first_spar_u
% 			l_interval_left_u(i) = sqrt((x_u(i)-x_u(i-1))^2 + (z_u(i)-z_u(i-1))^2);
% 		end
% 		perimeter_norm_left_u = sum(l_interval_left_u);
% 
% 		perimeter_norm_left = perimeter_norm_left_u + perimeter_norm_left_u;
% 
% 		%perimeter 4
% 		for i = index_last_spar_u:length(x_u_tot)
% 			l_interval_right_u(i) = sqrt((x_u_tot(i)-x_u_tot(i-1))^2 + (z_u_tot(i)-z_u_tot(i-1))^2);
% 		end
% 		perimeter_norm_right_u = sum(l_interval_right_u);
% 
% 		for i = index_last_spar_l:length(x_l_tot)
% 			l_interval_right_l(i) = sqrt((x_l_tot(i)-x_l_tot(i-1))^2 + (z_l_tot(i)-z_l_tot(i-1))^2);
% 		end
% 		perimeter_norm_right_l = sum(l_interval_right_l);
% 
% 		perimeter_norm_right = perimeter_norm_right_u + perimeter_norm_right_l;
% 
% 
% 
% 		%% Root
% 		perimeter_root = perimeter_norm_tot*c_root_w;
% 		spacing_root = perimeter_root/ns_wing;
% 
% 		x_wing_root_u = (x_u - X_centroid)*c_root_w;
% 		z_wing_root_u = (z_u - Z_centroid)*c_root_w;
% 
% 
% 		x_wing_root_l = (x_l - X_centroid)*c_root_w;
% 		z_wing_root_l = (z_l - Z_centroid)*c_root_w;
% 
% 
% 		%position of the stringers
% 		x_stringers_w(1) = x_wing_root_u(index_last_spar_u);
% 		z_stringers_w(1) = z_wing_root_u(index_last_spar_u);
% 		j = 2;
% 		lg = 0;
% 		for i = index_last_spar_u:-1:2
% 			lg = sqrt((x_wing_root_u(i)-x_wing_root_u(i-1))^2 + (z_wing_root_u(i)-z_wing_root_u(i-1))^2) + lg;
% 			if lg > spacing_root
% 				x_stringers_w(j) = x_wing_root_u(i);
% 				z_stringers_w(j) = z_wing_root_u(i);
% 				lg = 0;
% 				j = j +1;
% 			end
% 		end
% 
% 		for i = 2:index_last_spar_l
% 			lg = sqrt((x_wing_root_l(i)-x_wing_root_l(i-1))^2 + (z_wing_root_l(i)-z_wing_root_l(i-1))^2) + lg;
% 			if lg > spacing_root
% 				x_stringers_w(j) = x_wing_root_l(i);
% 				z_stringers_w(j) = z_wing_root_l(i);
% 				lg = 0;
% 				j = j +1;
% 			end
% 		end
% 
% 		x_stringers_w(ns_wing) =  x_wing_root_l(index_last_spar_l);
% 		z_stringers_w(ns_wing) =  z_wing_root_l(index_last_spar_l);
% 
% 
% 
% 		%% Tip
% 
% 		%washout
% 
% 		for i = 1:length(x_l_tot) % applying the washout
% 			x_l_tot_tip(i) = x_l_tot(i);
% 			z_l_tot_tip(i) = z_l_tot(i);
% 			d_l = norm([x_l_tot_tip(i), z_l_tot_tip(i)]);
% 			angle_l = atan(z_l_tot_tip(i)/x_l_tot_tip(i));
% 			if abs(x_l_tot_tip(i)) <= 1e-9 && abs(z_l_tot_tip(i)) <= 1e-9
% 				angle_l = 0;
% 			end
% 			x_l_tot_tip(i) = d_l*cos(angle_l + twist_rad);
% 			z_l_tot_tip(i) = d_l*sin(angle_l + twist_rad);
% 		end
% 
% 
% 
% 		for i = 1:length(x_u_tot) % applying the washout
% 			x_u_tot_tip(i) = x_u_tot(i);
% 			z_u_tot_tip(i) = z_u_tot(i);
% 			d_u = norm([x_u_tot_tip(i), z_u_tot_tip(i)]);
% 			angle_u = atan(z_u_tot_tip(i)/x_u_tot_tip(i));
% 			if abs(x_u_tot_tip(i)) <= 1e-9 && abs(z_u_tot_tip(i)) <= 1e-9
% 				angle_u = 0;
% 			end
% 			x_u_tot_tip(i) = d_u*cos(angle_u + twist_rad);
% 			z_u_tot_tip(i) = d_u*sin(angle_u + twist_rad);
% 		end
% 
% 
% 
% 		%spar positions
% 
% 		index_last_spar_u_tip = find(abs(x_u_tot_tip-x_pos_last_spar)<0.005);
% 		index_last_spar_l_tip = find(abs(x_l_tot_tip-x_pos_last_spar)<0.005);
% 
% 		index_first_spar_u_tip = find(abs(x_u_tot_tip-x_pos_first_spar)<0.005);
% 		index_first_spar_l_tip = find(abs(x_l_tot_tip-x_pos_first_spar)<0.005);
% 
% 		x_l_tip = x_l_tot_tip(1:index_last_spar_l_tip);
% 		z_l_tip = z_l_tot_tip(1:index_last_spar_l_tip);
% 		x_u_tip = x_u_tot_tip(1:index_last_spar_u_tip);
% 		z_u_tip = z_u_tot_tip(1:index_last_spar_u_tip);
% 
% 
% 
% 		perimeter_tip = perimeter_norm_tot*c_tip_w;
% 		spacing_tip = perimeter_tip/ns_wing;
% 
% 		x_wing_tip_u = x_u_tip*c_tip_w - X_centroid*c_root_w;%w.r.t the centroid of the root
% 		z_wing_tip_u = z_u_tip*c_tip_w - Z_centroid*c_root_w;%w.r.t the centroid of the root
% 
% 
% 		x_wing_tip_l = x_l_tip*c_tip_w - X_centroid*c_root_w;%w.r.t the centroid of the root
% 		z_wing_tip_l = z_l_tip*c_tip_w - Z_centroid*c_root_w;%w.r.t the centroid of the root
% 
% 
% 
% 		%position of the stringers
% 		x_stringers_w_tip(1) = x_wing_tip_u(index_last_spar_u_tip);
% 		z_stringers_w_tip(1) = z_wing_tip_u(index_last_spar_u_tip);
% 		j = 2;
% 		lg = 0;
% 		for i = index_last_spar_u_tip:-1:2
% 			lg = sqrt((x_wing_tip_u(i)-x_wing_tip_u(i-1))^2 + (z_wing_tip_u(i)-z_wing_tip_u(i-1))^2) + lg;
% 			if lg > spacing_tip
% 				x_stringers_w_tip(j) = x_wing_tip_u(i);
% 				z_stringers_w_tip(j) = z_wing_tip_u(i);
% 				lg = 0;
% 				j = j +1;
% 			end
% 		end
% 
% 		for i = 2:index_last_spar_l_tip
% 			lg = sqrt((x_wing_tip_l(i)-x_wing_tip_l(i-1))^2 + (z_wing_tip_l(i)-z_wing_tip_l(i-1))^2) + lg;
% 			if lg > spacing_tip
% 				x_stringers_w_tip(j) = x_wing_tip_l(i);
% 				z_stringers_w_tip(j) = z_wing_tip_l(i);
% 				lg = 0;
% 				j = j +1;
% 			end
% 		end
% 
% 		x_stringers_w_tip(ns_wing) =  x_wing_tip_l(index_last_spar_l_tip);
% 		z_stringers_w_tip(ns_wing) =  z_wing_tip_l(index_last_spar_l_tip);
% 
% 
% 		%% Stringers' cross section area
% 
% 		I_xx = sum(z_stringers_w.^2); % I_xx represent I_xx divided by B(the section)
% 		I_zz = sum(x_stringers_w.^2);
% 		I_xz = sum(x_stringers_w.*z_stringers_w); % =0 if x and/or z are axis of symetry
% 
% 		Bsigma_yy = zeros(ns_wing, length(V_points));
% 		maxBsigma_yy = zeros(1, length(V_points));
% 
% 		index_max_Mx_w = find(max(abs(Mx_w)));
% 		Mx_w_max = Mx_w(index_max_Mx_w);
% 
% 		index_max_Mz_w = find(max(abs(Mz_w)));
% 		Mz_w_max = Mz_w(index_max_Mz_w);
% 
% 		Bsigma_yy = ((-Mz_w_max*I_xx - Mx_w_max*I_xz)*x_stringers_w + (Mx_w_max*I_zz + Mz_w_max*I_xz)*z_stringers_w)/(I_xx*I_zz - I_xz^2); % eq 16.18 Megson [N]
% 		index_maxBsigma_yy = find(max(Bsigma_yy));
% 		maxBsigma_yy = Bsigma_yy(index_maxBsigma_yy); % [N]
% 
% 		B_min = maxBsigma_yy/sigma_max;
% 
% 		x_pos_first_spar = X_ac*c_root_w/c_bar_w - X_centroid*c_root_w;
% 		x_pos_last_spar = (1 - aileron_chord_ratio)*c_root_w;
% 
% 		index_first_spar = find(abs(x_stringers_w-x_pos_first_spar)<0.03);
% 
% 		l_spar_1 = sqrt((x_stringers_w(index_first_spar(1))-x_stringers_w(index_first_spar(2)))^2+ (z_stringers_w(index_first_spar(1))-z_stringers_w(index_first_spar(2)))^2);
% 		l_spar_2 = sqrt((x_stringers_w(1)-x_stringers_w(end))^2+ (z_stringers_w(1)-z_stringers_w(end))^2);
% 
% 		%% wall length
% 		l_12_curve = perimeter_norm_left*c_root_w;
% 		l_12_straight = l_spar_1;
% 		l_24 = perimeter_norm_l*c_root_w;
% 		l_43_straight = l_spar_2;
% 		l_43_curve = perimeter_norm_right*c_root_w;
% 		l_31 = perimeter_norm_u*c_root_w;
% 
% 
% 		%cell 1
% 
% 		Area1 = trapz(x_u_tot(1:index_first_spar_u),z_u_tot(1:index_first_spar_u)) - trapz(x_l_tot(1:index_first_spar_l),z_l_tot(1:index_first_spar_l));
% 		Pe1 = l_12_curve + l_12_straight;
% 		index_stringers_first_spar = find(abs(x_stringers_w -  x_wing_root_u(index_first_spar_u))<0.0025);
% 		X_stringers_1 = x_stringers_w(index_stringers_first_spar(1) : index_stringers_first_spar(2));
% 		Z_stringers_1 = z_stringers_w(index_stringers_first_spar(1) : index_stringers_first_spar(2));
% 
% 		%cell 2
% 
% 		Area2 = trapz(x_u_tot(index_first_spar_u:index_last_spar_u),z_u_tot(index_first_spar_u:index_last_spar_u)) - trapz(x_l_tot(index_first_spar_l:index_last_spar_l),z_l_tot(index_first_spar_l:index_last_spar_l));
% 		Pe2 = l_12_straight + l_24 + l_43_straight + l_31;
% 		% X_stringers_2 = [x_stringers_w(1 : index_stringers_first_spar(1))
% 		% x_stringers_w(index_stringers_first_spar(2):end)]; Z_stringers_2 =
% 		% [z_stringers_w(1 : index_stringers_first_spar(1))
% 		% z_stringers_w(index_stringers_first_spar(2):end)];
% 
% 		X_stringers_31 = x_stringers_w(1 : index_stringers_first_spar(1)-1);
% 		Z_stringers_31 = z_stringers_w(1 : index_stringers_first_spar(1)-1);
% 
% 		X_stringers_24 = x_stringers_w(index_stringers_first_spar(2)+1:end);
% 		Z_stringers_24 = z_stringers_w(index_stringers_first_spar(2)+1:end);
% 
% 		%cell 3
% 
% 		Area3 = trapz(x_u_tot(index_last_spar_u:end),z_u_tot(index_last_spar_u:end)) - trapz(x_l_tot(index_last_spar_l:end),z_l_tot(index_last_spar_l:end));
% 		Pe3 = l_43_curve + l_43_straight;
% 
% 		% % shear flow due to torque
% 		%
% 		% for i = 1:length(V_points)
% 		%     A = [[Area1, Area2];  %Momentum balance
% 		%         [Pe1/Area1+l_12_straight/Area2,
% 		%         -l_12_straight/Area1-Pe2/Area2]]; %Twist rate 1&2
% 		%     b = [My_w(i)/2, 0]'; q_M = A\b; q1_M(i) = q_M(1); % positive
% 		%     clockwise q2_M(i) = q_M(2); % positive clockwise
% 		% end
% 		%
% 		%
% 		%tapering
% 		delta_x = abs(x_stringers_w - x_stringers_w_tip);
% 		delta_y = span/2;
% 		delta_z = abs(z_stringers_w - z_stringers_w_tip);
% 
% 		T_x = zeros(1,length(V_points));
% 		T_z = zeros(1,length(V_points));
% 		for j = 1:ns_wing
% 			Px(j) = maxBsigma_yy*delta_x(j)/delta_y;
% 			Pz(j) = maxBsigma_yy*delta_z(j)/delta_y;
% 		end
% 		%
% 		%
% 		index_max_Tx_w = find(max(abs(Tx_w)));
% 		Tx_w_max = Tx_w(index_max_Tx_w);
% 
% 		index_max_Tz_w = find(max(abs(Tz_w)));
% 		Tz_w_max = Tz_w(index_max_Tz_w);
% 
% 		T_x = Tx_w_max - sum(Px);
% 		T_z = Tz_w_max - sum(Pz);
% 		%
% 		%
% 		% q_b1 = -(T_x*I_xx - T_z*I_xz)/(I_xx*I_zz -
% 		% I_xz^2)*B_min*sum(X_stringers_1) - (T_z*I_zz - T_x*I_xz)/(I_xx*I_zz -
% 		% I_xz^2)*B_min*sum(Z_stringers_1); q_b2 = -(T_x*I_xx -
% 		% T_z*I_xz)/(I_xx*I_zz - I_xz^2)*B_min*sum(X_stringers_2) - (T_z*I_zz -
% 		% T_x*I_xz)/(I_xx*I_zz - I_xz^2)*B_min*sum(Z_stringers_2);
% 
% 
% 		%torsion
% 		mu_ref = 25.9e9; %shear modulus of aluminium
% 
% 
% 		%% Open shear flow
% 		%circ 1-2
% 		q_0_12_c(1) = 0;
% 
% 		for i = 2:length(X_stringers_1)
% 			q_0_12_c(i) = q_0_12_c(i-1) - (T_x*I_xx - T_z*I_xz)/(I_xx*I_zz - I_xz^2)*B_min*X_stringers_1(i) - (T_z*I_zz - T_x*I_xz)/(I_xx*I_zz - I_xz^2)*B_min*Z_stringers_1(i);
% 		end
% 
% 		%3-1
% 		q_0_31(1) = 0;
% 
% 		for i = 2:length(X_stringers_31)
% 			q_0_31(i) = q_0_31(i-1) - (T_x*I_xx - T_z*I_xz)/(I_xx*I_zz - I_xz^2)*B_min*X_stringers_31(i) - (T_z*I_zz - T_x*I_xz)/(I_xx*I_zz - I_xz^2)*B_min*Z_stringers_31(i);
% 		end
% 
% 		%straight 1-2
% 		q_0_12 = q_0_31(end) - (T_x*I_xx - T_z*I_xz)/(I_xx*I_zz - I_xz^2)*B_min*X_stringers_1(1) - (T_z*I_zz - T_x*I_xz)/(I_xx*I_zz - I_xz^2)*B_min*Z_stringers_1(1);
% 
% 		%2-4
% 		q_0_24(1) = q_0_12_c(end) + q_0_12;
% 
% 		for i = 2:length(X_stringers_24)
% 			q_0_24(i) = q_0_24(i-1) - (T_x*I_xx - T_z*I_xz)/(I_xx*I_zz - I_xz^2)*B_min*X_stringers_24(i) - (T_z*I_zz - T_x*I_xz)/(I_xx*I_zz - I_xz^2)*B_min*Z_stringers_24(i);
% 		end
% 
% 
% 		%4-3 straight
% 		q_0_43 = q_0_24(end);
% 
% 
% 
% 		%% correction
% 		syms q_corr_1 q_corr_2 q_corr_3 dthetax_x_t
% 
% 		%--------------------------- Twist rate compatibility
% 		%---------------------------
% 		l_bar_1 = l_12_curve + l_12_straight;
% 		l_bar_2 = l_12_straight + l_24 + l_43_straight + l_31;
% 		l_bar_3 = l_43_curve + l_43_straight;
% 
% 		l_bar_12 = l_12_straight;
% 		l_bar_43 = l_43_straight;
% 
% 		%Cell1
% 		eqn1 = dthetax_x_t == 1/(2*Area1 * mu_ref) * (- q_corr_2*l_bar_12 + q_corr_1*l_bar_1  + sum(q_0_12_c * l_12_curve) - q_0_12*l_bar_12);
% 		%Cell2
% 		eqn2 = dthetax_x_t == 1/(2*Area2 * mu_ref) * (- q_corr_1*l_bar_12 - q_corr_3*l_bar_43 + q_corr_2*l_bar_2 + q_0_12*l_bar_12 + sum(q_0_24 .* l_24) + q_0_43*l_bar_43 + sum(q_0_31 .* l_31));
% 		%Cell3
% 		eqn3 = dthetax_x_t == 1/(2*Area3 * mu_ref) * (- q_corr_2*l_bar_43 + q_corr_3*l_bar_3 - q_0_43*l_bar_43);
% 
% 
% 
% 		%Swept area of each segment
% 		A_h = 0;
% 		t1_curve_12 = 0;
% 		for i =1:length(q_0_12_c)-1
% 			X1 = [X_stringers_1(i) Z_stringers_1(i)];
% 			X2 = [X_stringers_1(i+1) Z_stringers_1(i+1)];
% 			Xp = [X_centroid Z_centroid]*c_root_w;
% 			p = abs( (X2(1) - X1(1))*(X1(2) - Xp(2)) - (X1(1) - Xp(1))*(X2(2) - X1(2)) )/( (X2(1) - X1(1))^2 + (X2(2) - X1(2))^2 )^(1/2);
% 
% 			A_h = A_h + 0.5*p*( (X2(1) - X1(1))^2 + (X2(2) - X1(2))^2 )^(1/2);
% 			t1_curve_12 = t1_curve_12 + 2*A_h*q_0_12_c(i);
% 		end
% 
% 		A_h = 0;
% 		t1_13 = 0;
% 		for i =1:length(q_0_31)-1
% 			X1 = [X_stringers_31(i) Z_stringers_31(i)];
% 			X2 = [X_stringers_31(i+1) Z_stringers_31(i+1)];
% 			Xp = [X_centroid Z_centroid]*c_root_w;
% 			p = abs( (X2(1) - X1(1))*(X1(2) - Xp(2)) - (X1(1) - Xp(1))*(X2(2) - X1(2)) )/( (X2(1) - X1(1))^2 + (X2(2) - X1(2))^2 )^(1/2);
% 
% 			A_h = A_h + 0.5*p*( (X2(1) - X1(1))^2 + (X2(2) - X1(2))^2 )^(1/2);
% 			t1_13 = t1_13 + 2*A_h*q_0_31(i);
% 		end
% 
% 		A_h = 0;
% 		t1_24 = 0;
% 		for i =1:length(q_0_24)-1
% 			X1 = [X_stringers_24(i) Z_stringers_24(i)];
% 			X2 = [X_stringers_24(i+1) Z_stringers_24(i+1)];
% 			Xp = [X_centroid Z_centroid]*c_root_w;
% 			p = abs( (X2(1) - X1(1))*(X1(2) - Xp(2)) - (X1(1) - Xp(1))*(X2(2) - X1(2)) )/( (X2(1) - X1(1))^2 + (X2(2) - X1(2))^2 )^(1/2);
% 
% 			A_h = A_h + 0.5*p*( (X2(1) - X1(1))^2 + (X2(2) - X1(2))^2 )^(1/2);
% 			t1_24 = t1_24 + 2*A_h*q_0_24(i);
% 		end
% 
% 		t1_straight_12 = -q_0_12*l_bar_12*X_stringers_1(end);
% 
% 		t1_43 = q_0_43*l_bar_43*X_stringers_24(end);
% 
% 		t1 = t1_curve_12 + t1_straight_12 + t1_13 + t1_24 + t1_43;
% 		t2 = 2*(Area1*q_corr_1 + Area2*q_corr_2 + Area3*q_corr_3);
% 		t3 = sum(x_stringers_w.*Pz);
% 		t4 =  -  sum(z_stringers_w.*Px);
% 
% 		eqn4 = Mx_w_max == t1 + t2 + t3 + t4;
% 
% 
% 		[q_corr_1, q_corr_2, q_corr_3, dthetax_x_t]=vpasolve(eqn1, eqn2, eqn3, eqn4,q_corr_1, q_corr_2, q_corr_3, dthetax_x_t);
% 
% 		q_corr1=double(q_corr_1);
% 		q_corr2=double(q_corr_2);
% 		q_corr3=double(q_corr_3);
% 		dthetax_x_t=double(dthetax_x_t);
% 
% 
% 
% 		%% Total shear flow
% 
% 		q_12_c = q_0_12_c + q_corr1;
% 		q_12 = q_0_12 - q_corr1 + q_corr2;
% 		q_24 = q_0_24 + q_corr2;
% 		q_43 = q_0_43 + q_corr2 - q_corr3;
% 		q_43_c = q_corr3;
% 		q_31 = q_0_31 + q_corr2;
% 
% 
% 		%% Maximum shear stress
% 
% 		q_max = max( [ max(abs(q_12_c))  max(abs(q_12))  max(abs(q_24)) max(abs(q_43))  max(abs(q_43_c)) max(abs(q_31)) ]);
% 
% 		% Minimum thickness
% 		t_min_w = q_max/tau_max;
% 
	end
end