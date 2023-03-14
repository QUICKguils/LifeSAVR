% Last update: 27/02/2023

% TODO:
% - Check if there are comment about this graph in the intermediate report
%   correction.

function placard_diagram(opts)
% PLACARD_DIAGRAM  Altitude-velocity dependency.
%
% Parameters:
%	opts: char {'p'|'w'}, optional
%		'p' -> Enable plots creation.
%		'w' -> Write data in external file.
%
% Example:
%   placard_diagram('pw') % Create plots and write data to file.

% Se w default opts to 'p' (plot, but do not write).
if nargin == 0
	opts = 'p';
end

%% Data.

% Conversion factors.
ft2m = 0.3048;

% Quantities given or directly calculated from the statement.
M_cruise = 0.86;           % Mach number during cruise.
h = 30e3 * ft2m;           % Cruise altitude [m].
[~, ~, T_sl, ~] = ISA(0);  % Temperature at sea level [K].

%% Curves calculation.

% Sample of altitudes from sea level to stratospheric limit [m].
hsample = 0:100:11e3;

% Calculate air properties for all the sampled altitudes.
[rhos, ~, Ts, as] = arrayfun(@ISA, hsample);

% TAS for constant cruise mach [m/s].
tas_M = M_cruise * as .* sqrt(Ts/T_sl);

% Dynamic pressure at cruise altitude (Cst).
[rho_cruise, ~, T_cruise, a_cruise] = ISA(h);
V_cruise = M_cruise * a_cruise;
Cst = rho_cruise * (V_cruise .*sqrt(T_cruise/T_sl))^2/2;
% TAS such that the drag remains constant with altitude [m/s].
tas_D = sqrt(2*Cst./rhos);

% Plot the curves.
if contains(opts, 'p')
	hold on;
	plot(tas_M, hsample./1e3);
	plot(tas_D, hsample./1e3);
	title("Placard diagram.");
	xlabel("True airspeed [m/s]");
	ylabel("Altitude [km]");
	grid;
end

% Write data to file.
if contains(opts, 'w')
	% Open file for writing.
	fid = fopen("Plot/placard_diagram.txt", "w");

	% Write headers.
	fprintf(fid, "Altitude (km)\tTAS (m/s)\n");

	% Write data.
	fprintf(fid, "%f\t%f\n", [hsample/1e3, tas_D]);

	% Close file.
	fclose(fid);
end
end
