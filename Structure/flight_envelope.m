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

% Design cruise speed and design dive speed [m/s].
V_c = D.Propu.M_c * a;
V_d = D.Propu.M_d * a;

% Find the flight envelope points.
Man = maneuver_envelope();
% Gust = gust_envelope();

% Plot the curves.
if contains(opts, 'p')
	plot_envelope(Man.env.n, Man.env.EAS)
end

% % Write the data to external file.
% if contains(opts, 'w')
% 	write_ext();
% end

%% Maveuver envelope

	function Maneuver = maneuver_envelope()
		% MANEUVER_ENVELOPE  COmpute the maneuver envelope points.
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
		disp(n12);
		disp(V12);

		% Maximum upward load factor.
		V3 = V_d;
		n3 = C.n_up;

		% Weird interpolation.
		V4 = V_d;
		n4 = 0;

		% Maximum downward load factor.
		n5 = C.n_down;
		V5 = interp1([n4, -1], [V4, V_c], n5, "linear", "extrap");

		% Maximum negative lift.
		n6 = C.n_down;
		n61 = n1 + linspace(1, 0, 30).^2 * (n6 - n1);
		V61 = arrayfun(@v_max_lift, n61);

		% Conversion of true airspeed to equivalent airspeed (EAS).
		TAS2EAS = sqrt(rho/C.rho_sl);

		% Build Maneuver structure.
		% - Envelope: the envelope operating points, including stall
		%   lines. Useful for plotting.
		% - OP: contains only the envelope operating points. 
		Maneuver.env.n   = [n12, n3, n4, n5, n61];
		Maneuver.env.EAS = [V12, V3, V4, V5, V61] .* TAS2EAS;
		Maneuver.OP.n    = [n1,      n2,       n3, n4, n5, n6     ];
		Maneuver.OP.EAS  = [V12(1),  V12(end), V3, V4, V5, V61(1) ] .* TAS2EAS;
	end

%% Gust envelope

	function [n, EAS] = gust_envelope()
		% GUST_ENVELOPE

		% In the following of this section, we decided to label all the six
		% edges of the manever envelope with numeric indexes ranging from 1 to
		% 6. The numbering starts at flight condition (0 m/s; 0 g), and cycles
		% clockwise.

		% 14 CFR 23 recomandations (before 2017) for the gust values.
		gust_V_c = [50, 50, 25]    * C.ft2m;  % [m/s]
		gust_V_d = [25, 25, 12.5]  * C.ft2m;  % [m/s]
		gust_alt = [0, 20e3, 50e3] * C.ft2m;  % [m]
		% Interpolate these values.
		gust_c(alt) = @(alt) interp1(alt, gust_alt, gust_V_c, "linear");
		gust_d(alt) = @(alt) interp1(alt, gust_alt, gust_V_d, "linear");

		% Positive gust line at Vc.

		% Positive gust line at Vd.

		% Negative gust line at Vc.

		% Negative gust line at Vd.

		% Gather coordinate points.
		V = [V12, V3, V4, V5, V61];
		n = [n12, n3, n4, n5, n61];

		% We want to represent the equivalent airspeed (EAS).
		EAS = sqrt(rho/rho_sl) .* V;

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

	function plot_envelope(n, EAS) % TODO: Add gust
		% PLOT_ENVELOPE

		% Conversion in ISoU if desired.
		if ISoU
			EAS = EAS ./ C.kn2ms;
			h   = h   ./ C.ft2m;
		end

		% Instantiate a figure object.
		figure('WindowStyle', 'docked');

		% Plot EAS(n).
		TealDark = [000, 112, 127]/256;  % Taken from ULiège graphic chart.
			fill(EAS, n, TealDark, "FaceAlpha", 0.1);

		% Dress the plot.
		if ISoU
			title(strcat("Maneuver envelope (altitude of ", num2str(h), " ft)"));
			xlabel("Equivalent airspeed (kn)");
		else
			title(strcat("Maneuver envelope (altitude of ", num2str(h), " m)"));
			xlabel("Equivalent airspeed (m/s)");
		end
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