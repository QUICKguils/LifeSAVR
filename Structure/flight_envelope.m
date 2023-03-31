function flight_envelope(h, opts)
% FLIGHT_ENVELOPE  Build the flight envelope.
%
% Arguments:
%	h: double, optional.
%	  Altitude at which the flight envelope is desired, in meter.
%     Default is the cruise altitude.
%	opts: char {'p', 'w'}, optional.
%	  'p' -> Enable plots creation.
%	  'w' -> Write plotting data in external file.
%     Default is 'p'.

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
% - Generate the plot.
switch nargin
	case 0
		h  = C.h_cr;
		opts = 'p';
	case 1
		opts = 'p';
end

%% Main solve

% Function global variables.
[rho, ~, ~, a] = ISA(h);       % Air properties at desired altitude.
TAS2EAS = sqrt(rho/C.rho_sl);  % Conversion from true airspeed to equivalent airspeed.
V_c = D.Propu.M_c * a;         % Design cruise speed [m/s].
V_d = D.Propu.M_d * a;         % Design dive   speed [m/s].

% Build the flight envelope.
% In our case, the GE is entirely contained in the ME, so FE = ME.
ME = maneuver_envelope();
GE = gust_envelope();

% Retrieve relevant design speeds.
Speeds = design_EAS(ME, GE);

% Plot the curves.
if contains(opts, 'p')
	plot_envelope(ME, GE);
end

% % Write the data to external file.
% if contains(opts, 'w')
% 	write_ext();
% end

% Collect the relevant quantities to save, in a structure named `FE`.
FE.contour = ME.contour;
FE.CP      = ME.CP;
FE.speed   = Speeds;

% Save FE in data.mat, which lies in the root directory.
save(fullfile(file_dir, "../data.mat"), "FE", "-append");

%% Stall line points

	function TAS = tas_max_lift(n)
		% N_MAX_LIFT  True airspeed at max lift, for the given load factor.

		TAS = sqrt(2 * abs(n) * D.Plane.MTOW * C.g / (rho * D.Wing.surf * D.Wing.CL_max));
	end

%% Maveuver envelope

	function ME = maneuver_envelope()
		% MANEUVER_ENVELOPE  Build the maneuver envelope.

		% In this function, we decided to label all the six edges of the
		% maneuver envelope with numeric indexes ranging from 1 to 6. The
		% numbering starts at flight condition (0 m/s; 0 g), and cycles
		% clockwise.

		% Positive stall line.
		n_1    = 0;
		n_2    = C.n_up;
		n_12   = n_1 + linspace(0, 1, 30).^2 * (n_2 - n_1);  % Quadratically distributed.
		TAS_12 = arrayfun(@tas_max_lift, n_12);

		% Dive speed, maximum upward load factor.
		n_3   = C.n_up;
		TAS_3 = V_d;

		% Dive speed, free fall.
		n_4   = 0;
		TAS_4 = V_d;

		% Weird interpolation.
		n_5   = C.n_down;
		TAS_5 = interp1([n_4, -1], [TAS_4, V_c], n_5, "linear", "extrap");

		% Negative stall line.
		n_6    = C.n_down;
		n_61   = n_1 + linspace(1, 0, 30).^2 * (n_6 - n_1);  % Quadratically distributed.
		TAS_61 = arrayfun(@tas_max_lift, n_61);

		% Maneuver envelope structure.
		% - contour: sample of envelope contour points. Especially useful for plotting.
		% - CP:      envelope critical points, i.e. points 2 to 6.
		ME.contour = table(...
			[n_12,   n_3,   n_4,   n_5,   n_61  ]', ...
			[TAS_12, TAS_3, TAS_4, TAS_5, TAS_61]' .* TAS2EAS, ...
			'VariableNames', {'n', 'EAS'});
		ME.CP = table(...
			[n_2,         n_3,   n_4,   n_5,   n_6      ]', ...
			[TAS_12(end), TAS_3, TAS_4, TAS_5, TAS_61(1)]' .* TAS2EAS, ...
			'VariableNames', {'n', 'EAS'});
	end

%% Gust envelope

	function GE = gust_envelope()
		% GUST_ENVELOPE  Build the gust envelope.

		% 1. Gust lines.
		%
		% NOTE:
		% The computed gust lines should be more open that what is currently
		% obtained, leading to a gust envelope that is not entirely contained in
		% the maneuver envelope. However, calculations in this section has been
		% revised again and again, both with imperial and SI units. Reference
		% formulae can be found at:
		% www.ecfr.gov/on/2017-08-29/title-14/chapter-I/subchapter-C/part-23/subpart-C/

		% 14 CFR 23 recomandations (before 29-08-2017) for the gust values.
		Ue_c_alt = [50, 50,   25];    % [ft/s]
		Ue_d_alt = [25, 25,   12.5];  % [ft/s]
		gust_alt = [0,  20e3, 50e3];  % [ft]
		% Interpolate the gust speeds Ue for the desired altitude h.
		Ue_c = interp1(gust_alt, Ue_c_alt, h/C.ft2m, "linear");
		Ue_d = interp1(gust_alt, Ue_d_alt, h/C.ft2m, "linear");

 		% Airplaine weight ratio and gust alleviation factor.
		% FIX: we should use the CL_alpha of the plane, not the wing.
		mu = 2 * D.Plane.MTOW * C.g / (rho * D.Wing.CL_alpha * D.Wing.mac * C.g * D.Wing.surf);
		F = 0.88 * mu / (5.3 + mu);

		% Positive and negative gust lines at Vc and Vd.
		% EAS in expressed in [m/s].
		ng_c_plus_EAS  = @(EAS) 1 + F*C.rho_sl*D.Wing.surf*D.Wing.CL_alpha*(Ue_c*C.ft2m)*(EAS) / (2*D.Plane.MTOW*C.g);
		ng_c_minus_EAS = @(EAS) 1 - F*C.rho_sl*D.Wing.surf*D.Wing.CL_alpha*(Ue_c*C.ft2m)*(EAS) / (2*D.Plane.MTOW*C.g);
		ng_d_plus_EAS  = @(EAS) 1 + F*C.rho_sl*D.Wing.surf*D.Wing.CL_alpha*(Ue_d*C.ft2m)*(EAS) / (2*D.Plane.MTOW*C.g);
 		ng_d_minus_EAS = @(EAS) 1 - F*C.rho_sl*D.Wing.surf*D.Wing.CL_alpha*(Ue_d*C.ft2m)*(EAS) / (2*D.Plane.MTOW*C.g);
		% Convert the functions argument, to work with TAS.
		ng_c_plus_TAS  = @(TAS) ng_c_plus_EAS(TAS  * TAS2EAS);
		ng_c_minus_TAS = @(TAS) ng_c_minus_EAS(TAS * TAS2EAS);
		ng_d_plus_TAS  = @(TAS) ng_d_plus_EAS(TAS  * TAS2EAS);
 		ng_d_minus_TAS = @(TAS) ng_d_minus_EAS(TAS * TAS2EAS);

		% 2. Gust envelope.
		%
		% In this function, we decided to label all the six edges of the gust
		% envelope with numeric indexes ranging from 1 to 6. The numbering
		% starts at the minimum speed edge, and cycles clockwise.

		% Declare symbolic variable for speed, to solve gust lines and stall
		% lines intersections.
		TAS = sym('TAS');

		% Point 1: intersections of ng_c_minus with positive stall line.
		TAS_1 = double(vpasolve(tas_max_lift(ng_c_minus_TAS(TAS)) == TAS));
		n_1   = ng_c_minus_TAS(TAS_1);
		
		% Point 2: intersections of ng_c_plus with positive stall line.
		TAS_2 = double(vpasolve(tas_max_lift(ng_c_plus_TAS(TAS)) == TAS));
		n_2   = ng_c_plus_TAS(TAS_2);

		% Positive stall line.
		n_12    = n_1 + linspace(0, 1, 30).^2 * (n_2 - n_1);  % Quadratically distributed.
		TAS_12 = arrayfun(@tas_max_lift, n_12);

		% Point 3: positive V_c gust line at V_c.
		TAS_3 = V_c;
		n_3   = ng_c_plus_TAS(TAS_3);

		% Point 4: positive V_d gust line at V_d.
		TAS_4 = V_d;
		n_4   = ng_d_plus_TAS(TAS_4);

		% Point 5: negative V_d gust line at V_d.
		TAS_5 = V_d;
		n_5   = ng_d_minus_TAS(TAS_5);

		% Point 6: negative V_c gust line at V_c.
		TAS_6 = V_c;
		n_6   = ng_c_minus_TAS(TAS_6);

		% 3. Gust envelope structure.
		%
		% - contour: sample of envelope contour points. Especially useful for plotting.
		% - CP:      envelope critical points, i.e. points 1 to 6.
		% - line:    gust lines function handles.

		GE.contour = table(...
			[n_12,   n_3,   n_4,   n_5,   n_6,   n_1  ]', ...
			[TAS_12, TAS_3, TAS_4, TAS_5, TAS_6, TAS_1]' .* TAS2EAS, ...
			'VariableNames', {'n', 'EAS'});
		GE.CP = table(...
			[n_1,   n_2,   n_3,   n_4,   n_5,   n_6  ]', ...
			[TAS_1, TAS_2, TAS_3, TAS_4, TAS_5, TAS_6]' .* TAS2EAS, ...
			'VariableNames', {'n', 'EAS'});
		GE.line.ng_c_plus  = ng_c_plus_EAS;
		GE.line.ng_c_minus = ng_c_minus_EAS;
		GE.line.ng_d_plus  = ng_d_plus_EAS;
		GE.line.ng_d_minus = ng_d_minus_EAS;
	end

%% Design speeds

	function Speeds = design_EAS(ME, GE)
		% DESIGN_SPEEDS  Retrieve the relevant design EAS.
		%
		%   All the design speeds returned by this function have already been
		%   calculated. This function just nicely pack them in a dictionary.

		EAS_S = tas_max_lift(1) * TAS2EAS;  % Stall speed at n = 1 (cruise).
		EAS_A = ME.CP{1, 'EAS'};            % Stall speed at n = n_up.
		EAS_B = GE.CP{2, 'EAS'};            % Stall speed at n = ng_c_plus.
		EAS_C = V_c * TAS2EAS;              % Design cruise speed.
		EAS_D = V_d * TAS2EAS;              % Design dive speed.

		Speeds = dictionary(...
			["S",   "A",   "B",   "C",   "D"], ...
			[EAS_S, EAS_A, EAS_B, EAS_C, EAS_D]);
	end
%% Plot envelope

	function plot_envelope(ME, GE)
		% PLOT_ENVELOPE  Plot the ME, GE and gust lines.

		% Instantiate a figure object.
		figure('WindowStyle', 'docked');
		hold on;

		% Colors taken from ULiège graphic chart.
		TealDark   = [000, 112, 127] / 256;
		OrangeDark = [240, 127, 060] / 256;
		PurpleDark = [091, 037, 125] / 256;

		% Plot maneuver envelope.
		fill(ME.contour{:, 'EAS'} ./ C.kn2ms, ME.contour{:, 'n'}, TealDark, "FaceAlpha", 0.1);

		% Plot gust envelope.
		fill(GE.contour{:, 'EAS'} ./ C.kn2ms, GE.contour{:, 'n'}, PurpleDark, "FaceAlpha", 0.1);

		% Plot gust lines.
		v_lim = [0, 1.05 * V_d*TAS2EAS/C.kn2ms];
		fplot(@(EAS_kn) GE.line.ng_c_plus(EAS_kn*C.kn2ms),  v_lim, 'linestyle', '-.', 'Color', OrangeDark);
		fplot(@(EAS_kn) GE.line.ng_c_minus(EAS_kn*C.kn2ms), v_lim, 'linestyle', '-.', 'Color', OrangeDark);
		fplot(@(EAS_kn) GE.line.ng_d_plus(EAS_kn*C.kn2ms),  v_lim, 'linestyle', '-.', 'Color', OrangeDark);
		fplot(@(EAS_kn) GE.line.ng_d_minus(EAS_kn*C.kn2ms), v_lim, 'linestyle', '-.', 'Color', OrangeDark);

		% Dress the plot.
		title(strcat("Flight envelope (altitude of ", num2str(h/C.ft2m), " ft)"));
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