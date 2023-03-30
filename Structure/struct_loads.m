function structural_loads(AeroLoads)
% LOADS_STRUCTURAL  Structural loads.
%
% This function computes the structural loads exerting on the wings and
% fuselage cross sections, for all the critical points of the flight
% envelope.
%
% This function implements the methodology and formulae that can be
% found at:
% Aircraft Structures>lesson 6>slides 26 to 39.
%
% Argument:
%   AeroLoads: table
%	  Aerodynamic loads exerting on the wings, fuselage and tail, for
%	  all the critical points of the flight envelope.
%     This table is returned by Structure/loads_aero.m.

%% Imports

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Constants and aircraft data.
C = load(fullfile(file_dir, "../constants.mat"));
D = load(fullfile(file_dir, "../data.mat"));

% Utilities.
addpath(genpath(fullfile(file_dir, "../Utils")));

%% Main solve

% Compute MNT on the fuselage cross sections.
fuselage_loads();

% Compute MNT on the wings cross sections.
wings_loads();


%% Fuselage loads

	function fuselage_loads
		% FUSELAGE_LOADS  Loads on a fuselage cross section.

		L_fuselage = 8;
		L_engine = 1.74;

		x_rear_fuselage = CG_x(16) + 3/4*W.c_r;
		l_rear_fuselage = L_fuselage - x_rear_fuselage;
		W_rear_fuselage = l_rear_fuselage * Weight(1) * C.g;


		l_tail_fin = CG_x(15) - x_rear_fuselage;
		W_tail_fin = (Weight(14) + Weight(15))*C.g;

		l_engine = CG_x(20) - x_rear_fuselage;
		W_engine = Weight(20)*C.g;


		x_dd = L_fuselage - x_rear_fuselage;
		x_bb = (CG_x(20) - L_engine/2) - x_rear_fuselage;  %beginning of the engine
		x_cc = x_bb/2;

		D = [D_aa D_bb D_cc D_dd]; %[D_A D_B D_C D_D], D is constant over the entire rear fuselage length

		A_skin = pi*(max(D)+min(D))/2*l_rear_fuselage;

		W_maa = W_rear_fuselage * D_aa*pi/A_skin;
		W_mbb = W_rear_fuselage * D_bb*pi/A_skin;
		W_mcc = W_rear_fuselage * D_cc*pi/A_skin;
		W_mdd = W_rear_fuselage * D_dd*pi/A_skin;

		for i = 1:numel(V)
			SF_A(i) = n(i) * (W_mdd*l_rear_fuselage + (W_maa-W_mdd)*l_rear_fuselage/2 + W_tail_fin + W_engine);
			BM_A(i) = n(i)*cos(deg2rad(result_alpha(i))-0.0159)*(l_rear_fuselage^2/2*W_mdd + l_rear_fuselage^2/6*(W_maa-W_mdd) + l_tail_fin*W_tail_fin * L_engine*W_engine);
			% 0.0159 is the angle of incidence

			SF_C(i) = n(i)*(W_mdd*(l_rear_fuselage-x_cc) + (W_mcc-W_dd)*(l_rear_fus-x_cc)/2  + W_tail_fin + W_engine);
			BM_C(i) = n(i)*cos(deg2rad(result_alpha(i))-0.0159)*((l_rear_fuselage-x_cc)^2/2*W_mdd + (l_rear_fuselage-x_cc)^2/6*(W_mcc-W_mdd) + (l_tail_fin-x_cc)*W_tail_fin * (l_engine-x_cc)*W_engine);

			SF_B(i) = n(i)*(W_mdd*(l_rear_fuselage-x_bb) + (W_mcc-W_mdd)*(l_rear_fus-x_cc)/2 + W_tail_fin + W_engine);
			BM_B(i) = n(i)*cos(deg2rad(result_alpha(i))-0.0159)*((l_rear_fuselage-x_bb)^2/2*W_mdd + (l_rear_fuselage-x_bb)^2/6*(W_mbb-W_mdd)+ (l_tail_fin-x_bb)*W_tail_fin * (l_engine-x_bb)*W_engine);


			Ty_aa(i) = -F_fin(i);
			Tz_aa(i) = -(SF_A(i) - result_L_t(i))*cos(deg2rad(result_alpha(i))-0.0159);
			My_aa(i) = BM_A(i) - result_L_t(i)*l_tail_fin*cos(deg2rad(result_alpha(i))-0.0159);
			Mz_aa(i) = F_fin(i)*l_tail_fin;
			Mx_section_aa(i) = -F_fin(i)*abs(CG_z(14));


			Ty_cc(i) = -F_fin(i);
			Tz_cc(i) = -(SF_C(i) - result_L_t(i))*cos(deg2rad(result_alpha(i))-0.0159);
			My_cc(i) = BM_C(i) - result_L_t(i)*(l_tail_fin-x_cc)*cos(deg2rad(result_alpha(i))-0.0159);
			Mz_cc(i) = F_fin(i)*(l_tail_fin-x_cc);
			Mx_section_cc(i) = - F_fin(i)*abs(CG_z(14));


			Ty_bb(i) = -F_fin(i);
			Tz_bb(i) = -(SF_B(i) - result_L_t(i))*cos(deg2rad(result_alpha(i))-0.0159);
			My_bb(i) = BM_C(i) - result_L_t(i)*(l_tail_fin-x_bb)*cos(deg2rad(result_alpha(i))-0.0159);
			Mz_bb(i) = F_fin(i)*(l_tail_fin-x_bb);
			Mx_section_bb(i) = - F_fin(i)*abs(CG_z(14));
		end
	end

%% Wings loads

	function wings_loads
		% WINGS_LOADS  Loads on a wing cross section.

		for i = 1:numel(V)
			SF_A(i) = n(i) * (W_mdd*l_rear_fuslage + (W_maa-W_mdd)*l_rear_wing/2 + W_tail_fin + W_engine);
			BM_A(i) = n(i)*cos(deg2rad(result_alpha(i))-0.0159)*(l_rear_wing^2/2*W_mdd + l_rear_wing^2/6*(W_maa-W_mdd) + l_tail_fin*W_tail_fin * L_engine*W_engine);
			% 0.0159 is the angle of incidence

			SF_C(i) = n(i)*(W_mdd*(l_rear_wing-x_cc) + (W_mcc-W_dd)*(l_rear_fus-x_cc)/2  + W_tail_fin + W_engine);
			BM_C(i) = n(i)*cos(deg2rad(result_alpha(i))-0.0159)*((l_rear_wing-x_cc)^2/2*W_mdd + (l_rear_wing-x_cc)^2/6*(W_mcc-W_mdd) + (l_tail_fin-x_cc)*W_tail_fin * (l_engine-x_cc)*W_engine);

			SF_B(i) = n(i)*(W_mdd*(l_rear_wing-x_bb) + (W_mcc-W_mdd)*(l_rear_fus-x_cc)/2 + W_tail_fin + W_engine);
			BM_B(i) = n(i)*cos(deg2rad(result_alpha(i))-0.0159)*((l_rear_wing-x_bb)^2/2*W_mdd + (l_rear_wing-x_bb)^2/6*(W_mbb-W_mdd)+ (l_tail_fin-x_bb)*W_tail_fin * (l_engine-x_bb)*W_engine);


			Ty_aa(i) = -F_fin(i);
			Tz_aa(i) = -(SF_A(i) - result_L_t(i))*cos(deg2rad(result_alpha(i))-0.0159);
			My_aa(i) = BM_A(i) - result_L_t(i)*l_tail_fin*cos(deg2rad(result_alpha(i))-0.0159);
			Mz_aa(i) = F_fin(i)*l_tail_fin;
			Mx_section_aa(i) = -F_fin(i)*abs(CG_z(14));


			Ty_cc(i) = -F_fin(i);
			Tz_cc(i) = -(SF_C(i) - result_L_t(i))*cos(deg2rad(result_alpha(i))-0.0159);
			My_cc(i) = BM_C(i) - result_L_t(i)*(l_tail_fin-x_cc)*cos(deg2rad(result_alpha(i))-0.0159);
			Mz_cc(i) = F_fin(i)*(l_tail_fin-x_cc);
			Mx_section_cc(i) = - F_fin(i)*abs(CG_z(14));


			Ty_bb(i) = -F_fin(i);
			Tz_bb(i) = -(SF_B(i) - result_L_t(i))*cos(deg2rad(result_alpha(i))-0.0159);
			My_bb(i) = BM_C(i) - result_L_t(i)*(l_tail_fin-x_bb)*cos(deg2rad(result_alpha(i))-0.0159);
			Mz_bb(i) = F_fin(i)*(l_tail_fin-x_bb);
			Mx_section_bb(i) = - F_fin(i)*abs(CG_z(14));
		end
end