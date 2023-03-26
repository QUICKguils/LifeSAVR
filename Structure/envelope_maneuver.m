% TODO:
% - Reload data from propulsion.m
% - Something may be wrong with the chosen load factor limit. Recheck the KPP.
% - Recheck the gust envelope: it should not be contained in the maneuver
%   envelope !

function envelope_maneuver(h, opts)
% MANEUVER_ENVELOPE  Load_factor-velocity dependency.
%
% Parameter:
%	alt: double
%	  Altitude at which the maneuver envelope is desired, in meter.
%	opts: char {'p', 'w', 'i'}, optional
%		'p' -> Enable plots creation.
%		'w' -> Write data in external file.
%		'i' -> use the imperial system of units (abbr. ISoU).

close;

%% Data.

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% ISA function.
addpath(genpath(fullfile(file_dir, "Utils")));

% Conversion factors.
ft2m   = 0.3048;  % Foot -> meter.
kn2ms  = 0.5144;  % Knot -> meters per second.

% Option defaults:
% - Default altitude is the cruise altitude.
% - Yield the results in ISoU, and generate the plot.
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
[rho, ~, ~, a] = ISA(h);  % Air speed at desired altitude [m/s].

% Quantities retrieved from other parts.
CL_max = 1.493;  % Maximum lift coefficient.     - wing
S      = 7.37;   % Surface of the wing [m²].     - wing
MTOW   = 3353;   % Maximum takeoff weight [kg].  - weight

% Design cruise number and design dive Mach number.
% Lesson 4, slide 7 and 11.
M_C = 1.06 * M_cruise;
M_D = 1.07 * M_C;
% Design cruise speed and design dive speed [m/s].
V_C = M_C * a;
V_D = M_D * a;

% Determine if the imperial system of units is desired.
if contains(opts, 'i'); ISoU = true; else; ISoU = false; end

%% Curves calculations.
% In the following of this section, we decided to label all the six edges of the
% manever envelope with numeric indexes ranging from 1 to 6.
% The numbering starts at flight condition (0 m/s; 0 g), and cycles clockwise.

% Maximum positive lift.
n1 = 0;
n2 = n_up;
n12 = n1 + linspace(0, 1, 30).^2 * (n2 - n1);  % Quadratically distributed.
V12 = arrayfun(@v_max_lift, n12);
disp(n12);
disp(V12);

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
n61 = n1 + linspace(1, 0, 30).^2 * (n6 - n1);  % Quadratically distributed.
V61 = arrayfun(@v_max_lift, n61);

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

% Plot the curves.
if contains(opts, 'p')
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

% Write the data to external file.
if contains(opts, 'w')
	write_ext();
end

%% Local function definitions.

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

	% Print relevant velocities to the console output.
	% Stall TAS at cruise and at maximum upward load factor [m/s].
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