% TODO:
% - Assumptions on fuselage and wings specs are made a little bit
%   everywhere. Check that when more precise values are available.
% - Pass x_cs and y_cs as arguments, and define default values.
% - Wait for more precise CAD values.
% - Implement y option for wings_loads, if we have time.

function struct_loads(x_cs, y_cs)
% STRUCT_LOADS  Structural loads.
%
% This function computes the structural loads exerting on the wings and
% fuselage cross sections, for all the critical points of the flight
% envelope.
%
% This function implements the methodology and formulae that can be
% found at: Aircraft Structures>lesson 6>slides 26 to 39.
%
% Arguments:
%   x_cs: 1xn double, optional
%     X-coordinates of the desired fuselage cross sections [m].
%   y_cs: 1xn double, optional
%     Y-coordinates of the desired wing cross sections [m].
% Save:
%   FusLoads: table
%   WingLoads: table

%% Imports

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Constants and aircraft data.
C  = load(fullfile(file_dir, "../constants.mat"));
D  = load(fullfile(file_dir, "../data.mat"));
AL = D.AeroLoads;

% Utilities.
addpath(genpath(fullfile(file_dir, "../Utils")));

%% Options setting

% Option defaults (from the aircraft structures project statement):
% - cross section of the fuselage directly aft of the wing,
% - cross section of the wing at the root.
if ~nargin
	x_cs = D.Fus.x_start;
	y_cs = D.Fus.a;
end

%% Main solve

% Sanity checks on cross sections coordinates.
assert( ...
	all(x_cs >= D.Fus.x_start & x_cs <= D.Fus.x_end), ...
	"Desired fuselage cross section should lie in the rear fuselage.");
assert( ...
	all(y_cs >= D.Fus.a & y_cs <= D.Wing.span/2), ...  % TODO: not correct.
	"Desired wing cross section should lie in the wing.");

% The data structures to save consists of two table that contains the
% relevant MNT loads in the selected fuselage and wing cross sections,
% for al the CP.
FusLoads = table( ...
	'Size', [height(AL) * numel(x_cs), 8], ...
	'VariableNames', {'x',      'n',      'EAS',    'Ty',     'Tz',     'Mx',     'My',     'Mz'}, ...
	'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'});
WingLoads = table( ...
	'Size', [height(AL) * numel(y_cs), 8], ...
	'VariableNames', {'y',      'n',      'EAS',    'Tx',     'Tz',     'Mx',     'My',     'Mz'}, ...
	'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'});

% Fuselage pre-computations.
%
% Surface of the rear fuselage [m²].
S_rf = integral(@fuselage_perimeter, D.Fus.x_start, D.Fus.x_end, 'ArrayValued', true);
% Rear fuselage weight [N].
W_rf = D.Fus.mass * C.g;
% Aft fuselage linear weight [N/m].
w_end = fuselage_linweight(D.Fus.x_end);

% MNT in the selected fuselage cross sections.
nrow = 1;
for x = x_cs
	for idx = 1:height(AL)
		al = AL(idx, :);
		FusLoads(nrow, :) = fuselage_loads(al, x);
		nrow = nrow + 1;
	end
end

% MNT in the selected wings cross sections.
nrow = 1;
for y = y_cs
	for idx = 1:height(AL)
		al = AL(idx, :);
		WingLoads(nrow, :) = wings_loads(al);
		nrow = nrow + 1;
	end
end

% Save FusLoads and WingLoads in data.mat for the default cross sections
% values, otherwise simply display the results.
if ~nargin
	save(fullfile(file_dir, "../data.mat"), "FusLoads", "WingLoads", "-append");
else
	disp(FusLoads);
	disp(WingLoads);
end

%% Fuselage loads

	function loads = fuselage_loads(al, x)
		% FUSELAGE_LOADS  Loads on a fuselage cross section.

		% Angle of incidence of the plane [°].
		aoi = al.aoa - D.Wing.aoi;
		% Fuselage linear weight at cross section [N/m].
		w_x = fuselage_linweight(x);
		% Length from cross section to aft [m].
		L_rear = D.Fus.x_end - x;
		% Distance from HT to cross section [m].
		L_HT2cs = D.Comp{"Horizontal tail", "COG"}(1) - x;

		% Select the components that lies after the cross-section.
		comp = D.Comp(D.Comp.COG(:, 1) > x, :);
		% Retrieve their weights [N] and lever arms [m].
		larms   = comp.COG(:, 1) - x;
		weights = - comp.Mass * C.g;

		% Self weight of the fuselage is a trapezoïdal linear load.
		% We decompose it into a triangular linear load q1 and a constant
		% linear load q2.
		%
		% Equivalent loads Q1 and Q2 [N].
		Q1 = - w_end * L_rear;
		Q2 = -(w_x - w_end) * L_rear;
		% Equivalent lever arms [m].
		L1 = 1/3 * L_rear;
		L2 = 1/2 * L_rear;

		% Vertical dynamic equilibrium.
		F = al.n * (sum(weights) + Q1 + Q2);
		% Y-moments equilibrium.
		M = al.n * cosd(aoi) * (sum(weights.*larms) + Q1*L1 + Q2*L2);

		% Compute the relevant MNT loads.
		Ty = -al.F_fin;
		Tz = (F + al.P) * cosd(aoi);
		My = M + al.P * L_HT2cs * cosd(aoi);
		Mz = al.F_fin * L_HT2cs;
		Mx = -al.M_fus;  % See aero_loads.m: we don't take M_tail into account.

		% Return the computed loads.
		loads = table(x, al.n, al.EAS, Ty, Tz, Mx, My, Mz);
	end

	function perim = fuselage_perimeter(x)
		% FUSELAGE_PERIMETER  Perimeter of the fuselage along X-axis.
		%
		% Argument:
		%   x (double) -- X coordinate [m].
		% Return:
		%   perim (double) -- Perimeter [m].
		%
		% TODO: wait for precise fuselage geometrical specs.

		if x <= D.Fus.x_sectrans
			% Perimeter of an ellipe [m].
			perim = 2*pi * sqrt(D.Fus.a^2 + D.Fus.b^2 / 2);
		else
			% Perimeter of a circle [m].
			perim = pi * D.Fus.d;
		end
	end

	function w = fuselage_linweight(x)
		% FUSELAGE_LINWEIGHT  Linear weight of the fuselage along X-axis.
		%
		% Argument:
		%   x (double) -- X coordinate [m].
		% Return:
		%   w (double) -- Weight per unit length [N/m].
		%
		% TODO: wait for precise fuselage specs.

		% Perimeter of the rear fuselage at x [m].
		perim = fuselage_perimeter(x);

		w = W_rf * perim / S_rf;
	end

%% Wings loads
%
% WARN: Use lift, drag, etc. for ONE wing only.

	function loads = wings_loads(al)
		% WINGS_LOADS  Loads on a wing cross section.
		%
		% TODO:
		%  - This function is a draft, that only works for root cs.
		%  - We need more details about xyz coordinates of the components
		%    attached to the wing.

		% Angle of incidence of the plane [°].
		aoi = al.aoa - D.Wing.aoi;
		
		% Weight of one wing [N].
		W_wing = 0.5 * D.Wing.mass * C.g;
		% Coordinates of the wing COG [m].
		x_wing = D.Comp{"Wings", "COG"}(1);
		y_wing = D.Comp{"Wings", "COG"}(2);
		z_wing = D.Comp{"Wings", "COG"}(3);

		% Weight of the fuel stored in one wing [N].
		W_fuel = 0.5 * D.Comp{"Fuel Wing", "Mass"} * C.g;
		% Coordinates of the wing fuel COG [m].
		x_fuel = D.Comp{"Fuel Wing", "COG"}(1);
		y_fuel = D.Comp{"Fuel Wing", "COG"}(2);
		z_fuel = D.Comp{"Fuel Wing", "COG"}(3);

		% Compute the relevant MNT loads.
		% FIX: take into account the wing pitching moment.
		% FIX: choose reference origin from which to compute the moments.
		Tx = ( al.n*(W_wing + W_fuel) - al.L/2) * sind(aoi) + al.D_wing * cosd(aoi);
		Tz = (-al.n*(W_wing + W_fuel) + al.L/2) * cosd(aoi) + al.D_wing * sind(aoi);
		Mx = (-al.n*(W_wing*y_wing + W_fuel*y_fuel) - al.L/2*y_wing) * cosd(aoi);
		Mz = (-al.n*(W_wing*y_wing + W_fuel*y_fuel) - al.L/2*y_wing) * sind(aoi);
		My = (-al.n*(W_fuel*(x_fuel-x_wing))) * cosd(aoi);
		% TODO: better position in order to have correct Loads.
% 		My = (-al.n*(W_wing*x_wing + W_fuel*x_fuel) - al.L/2*x_wing + al.D_wing/2*x_wing) * cosd(aoi) ...
% 			+(-al.n*(W_wing*z_wing + W_fuel*z_fuel) - al.L/2*z_wing - al.D_wing/2*z_wing) * sind(aoi);

		% Return the computed loads.
		loads = table(y, al.n, al.EAS, Tx, Tz, Mx, My, Mz);
	end
end
