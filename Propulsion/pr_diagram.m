% Last update: 27/02/2023

% NOTE:
% - This file is based on slides 78-83 lesson 5 (Noels, coneptual design).
% - Actually, is exists an infinite number of payload-range diagram. We
%   should determine two of them: one for maximum range (by choosing the
%   optimal LD ratio), and the other for the maximum speed (ie for the
%   mach 0.86).
% - We should determine these graphs even for the fictious case of
%   replacing the payload by the fuel. (but not redesign the tank to
%   take into account this additional fuel capacity).

% TODO:
% - Wait for the optimal cruise to plot the max range diagram.
% - Recheck dry weight and fuel (volatile) weight.
% - reload data from propulsion.m. Engine has changed.
% - Check if there are comment about this graph in the intermediate report
%   correction.

function pr_diagram(condition, opts)
% PR_DIAGRAM  Plot of the payload-range diagram.
%
% Parameter:
%	condition: char {'s'|'r'}, optional
%		Specify the choice of flight condition.
%		's' -> PR diagram for maximum speed.
%		'r' -> PR diagram for maximum range.
%	opts: char {'p'|'w'}, optional
%		'p' -> Enable plots creation.
%		'w' -> Write data in external file.
%		'i' -> use the imperial system of units (abbr. ISoU).

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

close;

%% Data.

% Conversion factors.
ft2m  = 0.3048;
nmi2m = 1853.184;
lb2kg = 0.45359237;

% Physical constants.
g = 9.80665;  % Gravitational acceleration [m/s²].

% Quantities given or directly calculated from the statement.
t_loiter     = 5 * 3600;      % Minimum loiter time [s].
R_eg         = 2e3 * nmi2m;   % Minimum egress range [m].
R_ig         = 850 * nmi2m;   % Minimum ingress range [m].
M_cruise     = 0.86;          % Mach number during cruise.
M_loiter     = 0.7;           % Mach number during loiter.
h            = 30e3 * ft2m;   % Cruise altitude [m].
[~, ~, ~, a] = ISA(h);        % Sound speed [m/s].
V_loiter     = M_loiter * a;  % TAS at loiter [m/s].

% Quantities retrieved from other parts.
Cl         = 0.3;         % Coefficient of lift.         - wing
Cd         = 0.0225;      % Coefficient of drag.         - propu
SFC_cruise = 1.5188e-5;   % SFC at cruise [kg/(s*N)].    - propu
SFC_loiter = 1.4520e-5;   % SFC at loiter [kg/(s*N)].    - propu
W_MTO      = 3477;        % Maximum takeoff weight [kg]. - weight
W_payload  = 300 * lb2kg; % Maximum payload weight [kg]. - weight

%% Referral of all the weights.

% Total fuel weight in the tank [kg]. - propu
W_fuel_tot =  1315.3;

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
SFC_cruise = SFC_cruise * (1 + safety_bad_engine);
SFC_loiter = SFC_loiter * (1 + safety_bad_engine);

% Incemental weights [kg].
% In cruise, fuel weight for take-off is gone.
W_fuel = W_fuel_tot - W_fuel_takeoff;
% Fuel that stays in the plane (part of the dry weight).
W_fuel_dry = W_fuel_landing + W_fuel_trapped;
% Fuel that consumes during the flight (part of the volatile weight).
W_fuel_vol = W_fuel - W_fuel_dry;

% Total weights [kg].
% empty weight of the plane.
W_empty = W_MTO - W_fuel_tot - W_payload;
% Weight of the plane with payload, without the fuel.
W_nofuel = W_empty + W_payload;
% Dry weight of the plane: everithing except volatile fuel [kg].
W_dry = W_nofuel + W_fuel_dry;

%% Estimation of the maximum payload range.
% Obtained by summing the required ranges for the three main mission
% segments: ingress, loiter and egress.

% Equivalent loiter range [m].
R_loiter = t_loiter * V_loiter;
% Total mission range [m].
R_mission = R_ig + R_loiter + R_eg;
fprintf("Total range of the mission: %d km.\n", R_mission/1e3);

%% Payload-range for the desired fight condition.

[R, W_empty_R, W_nofuel_R, W_dry_R, W_tot_R] = payload_range(condition);

% Conversion in ISoU if desired.
if ISoU
	R          = R          ./ nmi2m;
	W_empty_R  = W_empty_R  ./ lb2kg;
	W_nofuel_R = W_nofuel_R ./ lb2kg;
	W_dry_R    = W_dry_R    ./ lb2kg;
	W_tot_R    = W_tot_R    ./ lb2kg;
end

disp(R(end));

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
	R = 1/SCF * eta_c * a/g * log((W_dry+W_vol)/W_dry);
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
		M  = M_cruise;
		SFC = SFC_cruise;
	elseif condition == 'r'
		M = M_loiter;
		SFC = SFC_loiter;
	end

	% Aerodynamic efficiency.
	eta_c = M * Cl/Cd;

	% Segment 1:
	% - Constant W_payload_1 = W_payload.
	% - Increase W_fuel_vol_1 from 0 to W_fuel_vol.
	R_1   = zeros(1, nsample);
	% Incremental weights.
	W_empty_1    = linspace(W_empty,    W_empty,    nsample);
	W_payload_1  = linspace(W_payload,  W_payload,  nsample);
	W_fuel_dry_1 = linspace(W_fuel_dry, W_fuel_dry, nsample);
	W_fuel_vol_1 = linspace(0,          W_fuel_vol, nsample);
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
	W_empty_2    = linspace(W_empty,    W_empty,    nsample);
	W_payload_2  = linspace(W_payload,  0,          nsample);
	W_fuel_dry_2 = linspace(W_fuel_dry, W_fuel_dry, nsample);
	W_fuel_vol_2 = linspace(W_fuel_vol, W_fuel_vol, nsample);
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
