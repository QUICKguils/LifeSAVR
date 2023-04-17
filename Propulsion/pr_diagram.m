% NOTE:
% - This file is based on slides 78-83 lesson 5 (Noels, coneptual design).
% - Actually, is exists an infinite number of payload-range diagram. We
%   should determine two of them: one for maximum range (by choosing the
%   optimal LD ratio), and the other for the maximum speed (ie for the
%   mach 0.86).

% TODO:
% - DO A TRADEOFF: make a bigger reservoir than necessary, to do a
%   payload/fuel tradeoff: show how the range evloves when we decide to
%   replace some payload by additional fuel.
% - Register stafety factors from propu in Propu structure.
% - Wait for the optimal cruise to plot the max range diagram.
% - Recheck dry weight and fuel (volatile) weight.

function pr_diagram(condition, opts)
% PR_DIAGRAM  Plot of the payload-range diagram.
%
% Arguments:
%	condition: char {'s'|'r'}, optional
%	  Specify the choice of flight condition.
%	  's' -> PR diagram for maximum speed.
%	  'r' -> PR diagram for maximum range.
%	opts: char {'p'|'w'}, optional
%	  'p' -> Enable plots creation.
%	  'w' -> Write data in external file.
%	  'i' -> use the imperial system of units (abbr. ISoU).

%% Imports

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Constants and aircraft data.
C = load(fullfile(file_dir, "../constants.mat"));
D = load(fullfile(file_dir, "../data.mat"));

%% Options setting

% Set argument defaults.
switch nargin
	case 0
		condition = 's';   % Maximum speed.
		opts      = 'pi';  % Enable plot, use ISoU.
	case 1
		opts = 'pi';
end

% Determine if the imperial system of units is desired.
if contains(opts, 'i'); ISoU = true; else; ISoU = false; end

%% Referral of all the weights.

% Total fuel weight in the tank [kg].
W_fuel_tot = D.Propu.W_fuel;

% Safety coefficients: account for extra fuel.
safety_takeoff      = 1e-2;  % Taxi and take-off.                    - propu
safety_landing      = 1e-2;  % Taxi and landing.                     - propu
safety_bad_engine   = 5e-2;  % Poorer-than-nominal fuel consumption. - propu
safety_trapped_fuel = 1e-2;  % Fuel trapped in the tank.             - propu
% Corresponding fuel weights [kg].
W_fuel_takeoff = W_fuel_tot * safety_takeoff;
W_fuel_landing = W_fuel_tot * safety_landing;
W_fuel_trapped = W_fuel_tot * safety_trapped_fuel;

% Fuel added for extra engine consumption can be seen as an increase in
% the engine SFC.
SFC_cruise = D.Propu.SFC_cruise * (1 + safety_bad_engine);
SFC_loiter = D.Propu.SFC_loiter * (1 + safety_bad_engine);

% Incemental weights [kg].
% In cruise, fuel weight for take-off is gone.
W_fuel = W_fuel_tot - W_fuel_takeoff;
% Fuel that stays in the plane (part of the dry weight).
W_fuel_dry = W_fuel_landing + W_fuel_trapped;
% Fuel that consumes during the flight (part of the volatile weight).
W_fuel_vol = W_fuel - W_fuel_dry;

% Total weights [kg].
% empty weight of the plane.
W_empty = D.Plane.MTOW - W_fuel_tot - C.W_payload;
% Weight of the plane with payload, without the fuel.
W_nofuel = W_empty + C.W_payload;
% Dry weight of the plane: everithing except volatile fuel [kg].
W_dry = W_nofuel + W_fuel_dry;

%% Estimation of the maximum payload range.
% Obtained by summing the required ranges for the three main mission
% segments: ingress, loiter and egress.

% Equivalent loiter range [m].
R_loiter = C.t_loiter_search * C.V_loiter;
% Total mission range [m].
R_mission = C.range_ingress + R_loiter + C.range_egress;

%% Payload-range for the desired fight condition.

[R, W_empty_R, W_nofuel_R, W_dry_R, W_tot_R] = payload_range(condition);

% Conversion in ISoU if desired.
if ISoU
	R          = R          ./ C.nmi2m;
	W_empty_R  = W_empty_R  ./ C.lb2kg;
	W_nofuel_R = W_nofuel_R ./ C.lb2kg;
	W_dry_R    = W_dry_R    ./ C.lb2kg;
	W_tot_R    = W_tot_R    ./ C.lb2kg;
end

% Plot the PR diagram.
if contains(opts, 'p')
	figure('WindowStyle', 'docked');
	hold on;
	plot(R, W_empty_R);
	plot(R, W_nofuel_R);
	plot(R, W_dry_R);
	plot(R, W_tot_R);
	if condition == 's'
		title('Payload-range diagram for maximum speed (M=0.86)');
	else
		title('Payload-range diagram for maximum range');
	end
	if ISoU
		xlabel('Range (nmi)');
		ylabel('Weight (lb)');
	else
		xlabel('Range (m)');
		ylabel('Weight (kg)');
	end
	grid;
	hold off;
end

% Write plotting data in an external file.
if contains(opts, 'w')
	write_ext();
end

%% Function definitions.

function R = range(SCF, eta_c, W_dry, W_vol)
	% RANGE  Estimate the range trough Breguet equation.
	%
	% Parameter:
	%	SFC_cruise: double
	%	  Specific fuel consumption of the engine at cruise [kg/(s*N)].
	%	eta_c: double
	%	  Aerodynamic efficiency.
	%	a, g: double
	%	  Speed of sound [m/s] and gravitational acceleration [m/s²].
	%	W_dry, W_vol: double
	%	  Dry weight and volatile weights of the plane [kg].
	% Return:
	%	R: double
	%	  Range of the plane [m].
	R = 1/SCF * eta_c * C.a_cr/C.g * log((W_dry+W_vol)/W_dry);
end

function [R, W_empty_R, W_nofuel_R, W_dry_R, W_tot_R] = payload_range(condition, nsample)
	% PAYLOAD_RANGE  Compute weights and associated ranges.
	%
	% Parameter:
	%	condition: char: {'s'|'r'}
	%		Flag that specify the choice of fight condition.
	%		's' -> PR diagram for maximum speed.
	%		'r' -> PR diagram for maximum range.
	%	nsample: double, optional
	%		Sample size for the weigths.
	% Return:
	%	R: double(1, 2*nsample)
	%		Ranges calculated for all the sampled weights.
	%	TOW, ZFC: double(1, 2*nsample)

	% Set default sample size to 20.
	if ~exist("nsample", "var")
		nsample = 20;
	end

	% Choose flight conditions:  max speed of max range.
	if condition == 's'
		M  = C.M_cr;
		SFC = SFC_cruise;
	elseif condition == 'r'
		M = C.M_loiter;
		SFC = SFC_loiter;
	end

	% Aerodynamic efficiency.
	eta_c = M * D.Wing.CL_cr/D.Propu.CD_plane;

	% Segment 1:
	% - Constant W_payload_1 = C.W_payload.
	% - Increase W_fuel_vol_1 from 0 to W_fuel_vol.
	R_1   = zeros(1, nsample);
	% Incremental weights.
	W_empty_1    = linspace(W_empty,     W_empty,     nsample);
	W_payload_1  = linspace(C.W_payload, C.W_payload, nsample);
	W_fuel_dry_1 = linspace(W_fuel_dry,  W_fuel_dry,  nsample);
	W_fuel_vol_1 = linspace(0,           W_fuel_vol,  nsample);
	% Total weights.
	W_nofuel_1 = W_empty_1  + W_payload_1;
	W_dry_1    = W_nofuel_1 + W_fuel_dry_1;
	W_tot_1    = W_dry_1    + W_fuel_vol_1;

	% Segment 2:
	% - Decrease TOW from MTOW to MTOW-MPW.
	% - Decrease PW from MPW to 0.
	% - Constant FW = MWF.
	R_2   = zeros(1, nsample);
	% Incremental weights.
	W_empty_2    = linspace(W_empty,     W_empty,    nsample);
	W_payload_2  = linspace(C.W_payload, 0,          nsample);
	W_fuel_dry_2 = linspace(W_fuel_dry,  W_fuel_dry, nsample);
	W_fuel_vol_2 = linspace(W_fuel_vol,  W_fuel_vol, nsample);
	% Total weights.
	W_nofuel_2 = W_empty_2  + W_payload_2;
	W_dry_2    = W_nofuel_2 + W_fuel_dry_2;
	W_tot_2    = W_dry_2    + W_fuel_vol_2;

	% Compute ranges.
	for i = 1:nsample
		R_1(i) = arrayfun(@range, SFC, eta_c, W_dry_1(i), W_fuel_vol_1(i));
		R_2(i) = arrayfun(@range, SFC, eta_c, W_dry_2(i), W_fuel_vol_2(i));
	end

	% Return values.
	R          = [R_1, R_2];
	W_empty_R  = [W_empty_1, W_empty_2];
	W_nofuel_R = [W_nofuel_1, W_nofuel_2];
	W_dry_R    = [W_dry_1, W_dry_2];
	W_tot_R    = [W_tot_1, W_tot_2];

end

function write_ext()
	% WRITE_EXT  Write plot data in an external file.

	% Specify the record file name.
	dirname = 'Plot/';
	% Name external file according to the units used.
	if ISoU; units = "-imperial"; else; units = "-metric"; end
	filename = strcat(... 
		dirname, ...
		'pr_diagram', ...
		'-condition_', condition, ...
		units);

	% Gather the data to store.
	pr_diagram_data = [R; W_empty_R; W_nofuel_R; W_dry_R; W_tot_R]';

	% Write in external file.
	writematrix(pr_diagram_data, filename);
end
end
