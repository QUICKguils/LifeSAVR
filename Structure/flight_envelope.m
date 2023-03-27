% TODO:
% - Reload data from propulsion.m
% - Something may be wrong with the chosen load factor limit. Recheck the KPP.
% - Recheck the gust envelope.

function flight_envelope(h, opts)
% FLIGHT_ENVELOPE
%
% TODO: write convention for flight points naming.
%
% Parameter:
%	alt: double
%	  Altitude at which the maneuver envelope is desired, in meter.
%	opts: char {'p', 'w', 'i'}, optional
%	  'p' -> Enable plots creation.
%	  'w' -> Write data in external file.
%	  'i' -> use the imperial system of units (abbr. ISoU).

%% Imports

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Constants and aircraft data.
C = load(fullfile(file_dir, "../constants.mat"));
D = load(fullfile(file_dir, "../data.mat"));

% Utilities.
addpath(genpath(fullfile(file_dir, "../Utils")));

%% Options setting

% Option defaults:
% - Default altitude is the cruise altitude.
% - Yield the results in ISoU, and generate the plot.
switch nargin
	case 0
		h  = C.h_cr;
		opts = 'pi';
	case 1
		opts = 'pi';
end

% Determine if the imperial system of units is desired.
if contains(opts, 'i'); ISoU = true; else; ISoU = false; end

%% Main solve

% Air properties at desired altitude.
[rho, ~, ~, a] = ISA(h);

% Conversion of true airspeed (TAS) to equivalent airspeed (EAS).
TAS2EAS = sqrt(rho/C.rho_sl);

% Design cruise speed and design dive speed [m/s].
V_c = D.Propu.M_c * a;
V_d = D.Propu.M_d * a;

% Find the flight envelope points.
Man = maneuver_envelope();
Gust = gust_envelope();
% Flight = flight_envelope(Man, Gust);

% Plot the curves.
if contains(opts, 'p')
	plot_envelope(Man, Gust);
end

% % Write the data to external file.
% if contains(opts, 'w')
% 	write_ext();
% end

%% Maveuver envelope

	function Maneuver = maneuver_envelope()
		% MANEUVER_ENVELOPE  Compute the maneuver envelope points.
		%
		% In this function, we decided to label all the six edges of the
		% maneuver envelope with numeric indexes ranging from 1 to 6. The
		% numbering starts at flight condition (0 m/s; 0 g), and cycles
		% clockwise.

		% Maximum positive lift.
		n1 = 0;
		n2 = C.n_up;
		% Quadratically distributed, to have smooth graph with a little
		% amount of points.
		n12 = n1 + linspace(0, 1, 30).^2 * (n2 - n1);
		V12 = arrayfun(@v_max_lift, n12);

		% Maximum upward load factor.
		V3 = V_d;
		n3 = C.n_up;

		% Dive speed, free fall.
		V4 = V_d;
		n4 = 0;

		% Weird interpolation.
		n5 = C.n_down;
		V5 = interp1([n4, -1], [V4, V_c], n5, "linear", "extrap");

		% Maximum negative lift.
		n6 = C.n_down;
		n61 = n1 + linspace(1, 0, 30).^2 * (n6 - n1);
		V61 = arrayfun(@v_max_lift, n61);

		% Build Maneuver structure.
		% - Envelope: the envelope operating points, including stall
		%   lines. Mainly useful for plotting.
		% - OP: contains only the envelope operating points. 
		Maneuver.env.n   = [n12, n3, n4, n5, n61];
		Maneuver.env.EAS = [V12, V3, V4, V5, V61] .* TAS2EAS;
		Maneuver.OP.n    = [n1,      n2,       n3, n4, n5, n6     ];
		Maneuver.OP.EAS  = [V12(1),  V12(end), V3, V4, V5, V61(1) ] .* TAS2EAS;
	end

%% Gust envelope

	function Gust = gust_envelope()
		% GUST_ENVELOPE  Build the gust lines functions.
		%
		% Reference values and formulae used to build the gust lines
		% can be found on:
		% www.ecfr.gov/on/2017-08-29/title-14/chapter-I/subchapter-C/part-23/subpart-C/

		% WARN: this function sporadically works with ISoU.
		W_ISoU   = D.Plane.MTOW * C.g / C.lbf2N;  % [lbf]
		rho_ISoU = rho / C.slug2kg * C.ft2m^3;    % [slug/ft³]
		mac_ISoU = D.Wing.mac / C.ft2m;           % [ft]
		g_ISoU   = C.g / C.ft2m;                  % [ft/s²]
		S_ISoU   = D.Wing.S / C.ft2m^2;           % [ft²]
		h_ISoU   = h / C.ft2m;                    % [ft]

		% 14 CFR 23 recomandations (before 29-08-2017) for the gust values.
		Ue_c_alt = [50, 50,   25];    % [ft/s]
		Ue_d_alt = [25, 25,   12.5];  % [ft/s]
		gust_alt = [0,  20e3, 50e3];  % [ft]
		% Interpolate the gust speeds Ue for the desired altitude h_ISoU.
		Ue_c = interp1(gust_alt, Ue_c_alt, h_ISoU, "linear");
		Ue_d = interp1(gust_alt, Ue_d_alt, h_ISoU, "linear");

		% Airplaine weight ratio and gust alleviation factor.
		% FIX: we should use the CL_alpha of the plane, not the wing.
		mu = 2 * W_ISoU / (rho_ISoU * D.Wing.CL_alpha * mac_ISoU * g_ISoU * S_ISoU);
		F = 0.88 * mu / (5.3 + mu);

		% Positive and negative gust lines at Vc and Vd.
		% WARN: unit: Ve is the airplane EAS, in knots.
		ng_c_plus  = @(Ve) 1 + F*D.Wing.CL_alpha*S_ISoU*Ue_c*Ve / (498 * W_ISoU);
		ng_c_minus = @(Ve) 1 - F*D.Wing.CL_alpha*S_ISoU*Ue_c*Ve / (498 * W_ISoU);
		ng_d_plus  = @(Ve) 1 + F*D.Wing.CL_alpha*S_ISoU*Ue_d*Ve / (498 * W_ISoU);
		ng_d_minus = @(Ve) 1 - F*D.Wing.CL_alpha*S_ISoU*Ue_d*Ve / (498 * W_ISoU);

		% The commented code below computes the gust lines with the SI
		% units. Unsurprizingly, wz obtain the same results.
		%
 		% % Airplaine weight ratio and gust alleviation factor.
		% % FIX: we should use the CL_alpha of the plane, not the wing.
		% mu = 2 * D.Plane.MTOW * C.g / (rho * D.Wing.CL_alpha * D.Wing.mac * C.g * D.Wing.S);
		% F = 0.88 * mu / (5.3 + mu);
		% % Positive and negative gust lines at Vc and Vd.
		% ng_c_plus  = @(Ve) 1 + F*C.rho_sl*D.Wing.S*D.Wing.CL_alpha*(Ue_c*C.ft2m)*(Ve*C.kn2ms) / (2*D.Plane.MTOW*C.g);
		% ng_c_minus = @(Ve) 1 - F*C.rho_sl*D.Wing.S*D.Wing.CL_alpha*(Ue_c*C.ft2m)*(Ve*C.kn2ms) / (2*D.Plane.MTOW*C.g);
		% ng_d_plus  = @(Ve) 1 + F*C.rho_sl*D.Wing.S*D.Wing.CL_alpha*(Ue_d*C.ft2m)*(Ve*C.kn2ms) / (2*D.Plane.MTOW*C.g);
 		% ng_d_minus = @(Ve) 1 - F*C.rho_sl*D.Wing.S*D.Wing.CL_alpha*(Ue_d*C.ft2m)*(Ve*C.kn2ms) / (2*D.Plane.MTOW*C.g);

		% Build Gust structure.
		Gust.ng_c_plus  = ng_c_plus;
		Gust.ng_c_minus = ng_c_minus;
		Gust.ng_d_plus  = ng_d_plus;
		Gust.ng_d_minus = ng_d_minus;
	end

%% Flight envelope

	function Flight = flight_envelope(Man, Gust)
		% FLIGHT_ENVELOPE  Compute the flight envelope.

		% Vs: intersection of stall line with n=1. 

		% Va: intersection of stall line with n=C.n_up
		% Vb: intersection of stall line with ng_c_plus

	end

%% Stall line points

	function TAS = v_max_lift(n)
		% N_MAX_LIFT  True airspeed at max lift, for the given load factor.
		%
		% Parameter:
		%	n: double
		%		Load factor.
		TAS = sqrt(2 * abs(n) * D.Plane.MTOW * C.g / (rho * D.Wing.S * D.Wing.CL_max));
	end

%% Plot envelope

	function plot_envelope(Man, Gust)
		% PLOT_ENVELOPE

		% Conversion in ISoU.
		Man.env.EAS = Man.env.EAS   ./ C.kn2ms;
		h_ISoU      = h             ./ C.ft2m;
		V_c_KEAS    = V_c * TAS2EAS ./ C.kn2ms;
		V_d_KEAS    = V_d * TAS2EAS ./ C.kn2ms;

		% Instantiate a figure object.
		figure('WindowStyle', 'docked');
		hold on;

		% Plot maneuver envelope.
		TealDark = [000, 112, 127]/256;  % Taken from ULiège graphic chart.
		fill(Man.env.EAS, Man.env.n, TealDark, "FaceAlpha", 0.1);

		% Plot gust lines.
		fplot(Gust.ng_c_plus,  [0, 1.05*V_d_KEAS], 'linestyle', '-.', 'Color', [0.9290 0.6940 0.1250]);
		fplot(Gust.ng_c_minus, [0, 1.05*V_d_KEAS], 'linestyle', '-.', 'Color', [0.9290 0.6940 0.1250]);
		fplot(Gust.ng_d_plus,  [0, 1.05*V_d_KEAS], 'linestyle', '-.', 'Color', [0.9290 0.6940 0.1250]);
		fplot(Gust.ng_d_minus, [0, 1.05*V_d_KEAS], 'linestyle', '-.', 'Color', [0.9290 0.6940 0.1250]);

		% Plot flight envelope.

		% Dress the plot.
		title(strcat("Flight envelope (altitude of ", num2str(h_ISoU), " ft)"));
		xlabel("Equivalent airspeed (kn)");
		ylabel("Load factor");
		grid;
		axis padded;
	end

%% Write plot to extenal file

% 	function write_ext()
% 		% Name external file according to the units used.
% 		if ISoU; units = "imperial"; else; units = "metric"; end
% 		filename = strcat("Plot/maneuver_envelope-", units, "-alt_", num2str(h), ".dat");
% 
% 		% Write data to file.
% 		writematrix([EAS; n;]', filename);
% 
% 		% Conversion in ISoU if desired.
% 		if ISoU
% 			EAS_st = EAS_st ./ kn2ms;
% 			EAS_A  = EAS_A  ./ kn2ms;
% 			EAS_C  = EAS_C  ./ kn2ms;
% 			EAS_D  = EAS_D  ./ kn2ms;
% 		end
% 		if ISoU; units = "kn"; else; units = "m/s"; end
% 		fprintf("Stall speed at cruise: %f " + units + "\n", EAS_st);
% 		fprintf("Maximum lift speed:    %f " + units + "\n", EAS_A);
% 		fprintf("Design cruise speed:   %f " + units + "\n", EAS_C);
% 		fprintf("Design dive speed:     %f " + units + "\n", EAS_D);
% 	end

end