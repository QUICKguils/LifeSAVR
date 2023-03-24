% Last update: 27/02/2023

% TODO:
% - Reload data from propulsion.m
% - Something may be wrong with the chosen load factor limit. Recheck the KPP.
% - Recheck the gust envelope: it should not be contained in the maneuver envelope !

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

close;

%% Data.

% Conversion factors.
ft2m   = 0.3048;  % Foot -> meter.
kn2ms  = 0.5144;  % Knot -> meters per second.

% Option defaults: - Default altitude is the cruise altitude. - Yield
% the results in ISoU, and generate the plot.
switch nargin
	case 0
		h  = 30e3 * ft2m;   % Cruise altitude [m].
		opts = 'pi';
	case 1
		opts = 'pi';
end

% Physical constants and directly computable quantities.
[rho_sl, ~, ~, ~] = ISA(0);  % Air density at sea level [kg/m³].
g = 9.80665;                 % Gravitational acceleration [m/s²].

% Quantities given or directly calculated from the statement.
n_up           = 3;       % KPP 18 - Upward flight limit load factor.
n_dw           = -1.5;    % KPP 19 - Downward flight limit load factor.
M_cruise       = 0.86;    % KPP 03 - Mach number during cruise.
[rho, ~, ~, a] = ISA(h);  % Air speed at cruise altitude [m/s].

% Quantities retrieved from other parts.
CL_max = 1.493;  % Maximum lift coefficient.     - wing
S      = 7.37;   % Surface of the wing [m²].     - wing
MTOW   = 3353;   % Maximum takeoff weight [kg].  - weight

% Design cruise number and design dive Mach number. Lesson 4, slide 7
% and 11.
M_C = 1.06 * M_cruise;
M_D = 1.07 * M_C;
% Design cruise speed and design dive speed [m/s].
V_C = M_C * a;
V_D = M_D * a;

% Determine if the imperial system of units is desired.
if contains(opts, 'i'); ISoU = true; else; ISoU = false; end

%% Main solve.

% Find the flight envelope points.
a = maneuver_envelope();
b = gust_envelope();

% Plot the curves.
if contains(opts, 'p')
	plot_envelope()
end

% Write the data to external file.
if contains(opts, 'w')
	write_ext();
end

%% Maveuver envelope.

	function maneuver_envelope()
		% MANEUVER_ENVELOPE
		%
		% In the following of this section, we decided to label all the six
		% edges of the manever envelope with numeric indexes ranging from 1 to
		% 6. The numbering starts at flight condition (0 m/s; 0 g), and cycles
		% clockwise.

		% Maximum positive lift.
		n0 = 0;
		nA = n_up;
		% Quadratically distributed, to have smooth graph with a little
		% amount of points.
		n0A = n1 + linspace(0, 1, 30).^2 * (nA - n0);
		V0A = arrayfun(@v_max_lift, n0A);

		% Maximum upward load factor.
		V3 = V_D;
		n3 = n_up;

		% Weird interpolation.
		V4 = V_D;
		n4 = 0;

		% Maximum downward load factor.
		n5 = n_dw;
		V5 = interp1([n4, -1], [V4, V_C], n5, "linear", "extrap");

		% Maximum negative lift.
		n6 = n_dw;
		n61 = n1 + linspace(1, 0, 30).^2 * (n6 - n1);
		V61 = arrayfun(@v_max_lift, n61);

		% Gather coordinate points.
		V = [V0A, V3, V4, V5, V61];
		n = [n0A, n3, n4, n5, n61];

		% We want to represent the equivalent airspeed (EAS).
		EAS = sqrt(rho/rho_sl) .* V;

		% Conversion in ISoU if desired.
		if ISoU
			EAS = EAS ./ kn2ms;
			h = h ./ ft2m;
		end
	end

%% Gust envelope.

	function gust_envelope()
		% GUST_ENVELOPE

		% In the following of this section, we decided to label all the six
		% edges of the manever envelope with numeric indexes ranging from 1 to
		% 6. The numbering starts at flight condition (0 m/s; 0 g), and cycles
		% clockwise.

		% 14 CFR 23 recomandations (before 2017) for the gust values
		gust_Vc = [50, 50, 25] * ft2m;    % [m/s]
		gust_Vd = [25, 25, 12.5] * ft2m;  % [m/s]
		gust_h = [0, 20e3, 50e3] * ft2m;  % [m]
		% Interpolate these values.

		% Positive gust line at Vc.

		% Positive gust line at Vd.

		% Negative gust line at Vc.

		% Negative gust line at Vd.

		% Gather coordinate points.
		V = [V12, V3, V4, V5, V61];
		n = [n12, n3, n4, n5, n61];

		% We want to represent the equivalent airspeed (EAS).
		EAS = sqrt(rho/rho_sl) .* V;

		% Conversion in ISoU if desired.
		if ISoU
			EAS = EAS ./ kn2ms;
			h = h ./ ft2m;
		end
	end

%% Local function definitions.

	function plot_envelope()
		% PLOT_ENVELOPE

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

	function TAS = v_max_lift(n)
		% N_MAX_LIFT  True airspeed at max lift, for the given load factor.
		%
		% Parameter:
		%	n: double
		%		Load factor.
		TAS = sqrt(2 * abs(n) * MTOW * g / (rho * S * CL_max));
	end

	function write_ext()
		% Name external file according to the units used.
		if ISoU; units = "imperial"; else; units = "metric"; end
		filename = strcat("Plot/maneuver_envelope-", units, "-alt_", num2str(h), ".dat");

		% Write data to file.
		writematrix([EAS; n;]', filename);

		% Print relevant velocities to the console output. Stall TAS at cruise
		% and at maximum upward load factor [m/s].
		V_st = v_max_lift(1);
		V_A = V12(end);
		% Conversion to EAS [m/s].
		EAS_st = sqrt(rho/rho_sl) .* V_st;
		EAS_A  = sqrt(rho/rho_sl) .* V_A;
		EAS_C  = sqrt(rho/rho_sl) .* V_C;
		EAS_D  = sqrt(rho/rho_sl) .* V_D;
		% Conversion in ISoU if desired.
		if ISoU
			EAS_st = EAS_st ./ kn2ms;
			EAS_A  = EAS_A  ./ kn2ms;
			EAS_C  = EAS_C  ./ kn2ms;
			EAS_D  = EAS_D  ./ kn2ms;
		end
		if ISoU; units = "kn"; else; units = "m/s"; end
		fprintf("Stall speed at cruise: %f " + units + "\n", EAS_st);
		fprintf("Maximum lift speed:    %f " + units + "\n", EAS_A);
		fprintf("Design cruise speed:   %f " + units + "\n", EAS_C);
		fprintf("Design dive speed:     %f " + units + "\n", EAS_D);
	end
end
