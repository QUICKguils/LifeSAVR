% PROPULSION  Calculations related to the propulsion.
%   Last update: 13/03/2023

% PROCEDURE:
% - Get an estimation of the thrust needed for the dash speed.
% - Translate this thrust to the equivalent thrust at sea level.
% - Take a safety of around 20% more. We get the installed thrust.
% - Translate this installed thrust to the uninstalled thrust. It
%   corresponds to the thrust measured in lab conditions by the
%   manufacturers.
% - Select a few (3-4) engines that match this uninstalled thrust.
% - Find the SFC specified by the manufacturers.
% - Translate this SFC to the cruise conditions.
% - We can now calculate the fuel weight, given the mission duration (no
%   need to take into account the takeoff/landing at this stage, or
%   maybe add someting like 2%).
% ----------
% - Next step, calculate the total weight, by summing all the parts
%   (wings, tails, landing gears, engine, ...).
% - And ONLY at this moment, we can see if the thrust will be sufficient
%   for the takeoff (normally it will be sufficient, as no climb rate is
%   specified, and the runway length is decent).
% - And now we are in an iterative loop: see if the total weight when
%   summing all the parts are correct (respect the total weight limit).
%   Else, redo the design loop of the aircraft. For this part (engine)
%   just put a different engine.

% TODO:
% - use sweep angle at LE in the estimation of e.
% - Check the 20% and the 8% factor: 20% seems too high. Provide refs.
% - Think more about the air duct. What size, how to beneficiate from
%   the BL of the nose/fuselage. See lesson 5 (Koen).
% - Carry out a small tradeoff study for the engine: BRP and SFC vs
%   engine diameter (fuselage drag).

clear;

%% Data.

% Conversion factors.
ft2m  = 0.3048;
nmi2m = 1853.184;
lb2kg = 0.45359237;
lbf2N = 4.4482216152605;

% Physical constants.
g = 9.80665;  % Gravitational acceleration [m/s²].

% Quantities given or directly calculated from the statement.
t_loiter1      = 5 * 3600;      % Minimum loiter 1 time [s].
t_loiter2      = 45 * 60;       % Minimum loiter 2 time [s].
eg_range       = 2e3 * nmi2m;   % Minimum egress range [m].
ig_range       = 850 * nmi2m;   % Minimum ingress range [m].
M_cruise       = 0.86;          % Mach number during cruise.
M_loiter       = 0.7;           % Mach number during loiter.
M_dash         = 0.9;           % Mach number at dash speed.
h              = 30e3 * ft2m;   % Cruise altitude [m].
[rho, p, T, a] = ISA(h);        % Fluid properties at altitude h. See ISA.m
V_cruise       = M_cruise * a;  % TAS at cruise [m/s].
V_loiter       = M_loiter * a;  % TAS at loiter [m/s].
V_dash         = M_dash   * a;  % TAS at dash   [m/s]
W_payload      = 300 * lb2kg;   % Payload weight [kg].

% Quantities retrieved from other parts.
Cl     = 0.276;  % Coefficient of lift.               - wing
AR     = 8;      % Aspect ratio.                      - wing
Lambda = 24;     % Quarter chord sweep angle [°].     - wing
e      = 0.612;  % Oswald's efficiency factor.        - wing
S      = 7.36;   % Surface of the wing [m²].          - wing
Cd_0   = 0.0172; % Lift-independent drag coefficient. - wing
MTOW   = 3367;   % Initial guess of MTOW [kg].        - aircraft

% Selected engine data.
% Thrust, SFC and BPR data are gathered by looking at either the
% manufacturer's datasheets or the engine certificate (FAA or EASA).
% G is the "gas generator function". Correlation between G and BPR can
% be found in:
% Teams>Propulsion>Files>Sources>Bartel-SFC_Calculations.pdf (page 5)
engine_table = table( ...
	'Size',          [0, 6], ...
	'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double'}, ...
	'VariableNames', {'Brand',  'Thrust', 'SFC',    'BPR',    'G',      'n'});

% Register selected engines.
% Models cited below belongs to engine series. Look at other models in
% the series if more thrust is needed.
engine_table("PW535E1",     :) = {"Pratt & Whitney",        15.47e3, 1.25e-5, 2.55, 0.78, 0.5};   % F = 1396 kg,  E = 317 kg, D = 1.08 m
engine_table("PW306B",      :) = {"Pratt & Whitney",        23.24e3, 1.16e-5, 4.5,  0.91, 0.8};   % F = 1442 kg,  E = 450 kg, D = 0.97 m, a bit more than 2 m long
engine_table("JT15D-5C",    :) = {"Pratt & Whitney",        14.18e3, 1.59e-5, 2,    0.75, 0.5};   % F = 1769 kg,  E = 302 kg, D = 0.69 m
engine_table("FJ44-4A",     :) = {"Williams International", 16.08e3, 1.31e-5, 3.4,  0.84, 0.75};  % F = 1688 kg,  E = 304 kg, D = 0.82 m
engine_table("TFE731-20BR", :) = {"Honeywell Aerospace",    16.2e3,  1.25e-5, 3.1,  0.82, 0.75};  % F = 1612,     E = 405 kg, D = 1 m (0.85 for the fan)
engine_table("CFE738",      :) = {"CFE",                    26.3e3,  1.05e-5, 5.3,  0.98, 0.8};   % F = 1398 kg,  E = 551,    D = 0.9 m >
engine_table("CFE700",      :) = {"CFE",                    20.24e3, 1.85e-5, 2,    0.75, 0.5};   % F = 20541 kg, E = 333,    D = 0.84 m
engine_table("HF120",       :) = {"GE Honda",               9.1e3,   1.99e-5, 2.9,  0.8,  0.5};

% Select one engine.
engine = table2struct(engine_table("FJ44-4A", :));

%% Estimation of the thrust needed for dash speed.
% Slides 27 and following, lesson 2.

% Estimation of e is now calculated in the wing part.
% % Estimation of the Oswald's efficiency factor.
% % Teams>General>Books>Raymer (page 444)
% e_straight = 1.78 * (1 - 0.045*AR^0.68)                     - 0.64;
% e_swept    = 4.61 * (1 - 0.045*AR^0.68) * cosd(Lambda)^0.15 - 3.1;
% % Linear interpolation for a sweep angle of Lambda.
% e = e_straight + (e_swept - e_straight)/30 * Lambda;

% Estimation of the drag coefficient.
Cd = Cd_0 + Cl^2 / (e * pi * AR);

% Thrust at cruise must equal the drag [N].
F_cruise = 0.5 * rho * V_cruise^2 * S * Cd;

%% Cruise thrust to uninstalled thrust.
% The uninstalled thrust corresponds to the thrust measured in lab
% conditions by the manufacturers.
% WARN: conditions assumed below can depend on the manufacturer!

% Lab flow conditions: assume sea level.
h_sl = 0;                               % Altitude (sea level) [m].
[rho_sl, p_sl, T_sl, a_sl] = ISA(h_sl); % Fluid properties at sea level. See ISA.m

% Thrust at sea level [N].
% Actually there is a causal dilemma: we want to estimate the equivalent
% uninstalled thrust to select a matching engine. But at the same time,
% we need to know which engine is used in order to evaluate this
% equivalent uninstalled thrust. More precisely, the equivalent thrust
% depend on the engine BPR.
% -> We can have a first guesstimate of 0.8 for the value of G. This
% corresponds to a quite small BPR engine.
% Teams>Propulsion>Files>Sources>Bartel-SFC_Calculations.pdf (page 5)
BPR = engine.BPR;
G = engine.G;
A = -0.4327 * (p/p_sl)^2 + 1.3855 * p/p_sl + 0.0472;
X =  0.1377 * (p/p_sl)^2 - 0.4374 * p/p_sl + 1.3003;
Z =  0.9106 * (p/p_sl)^2 - 1.7736 * p/p_sl + 1.8697;
F_sl = F_cruise / ( ...
	  A ...
	- Z * p/p_sl * (0.377*(1+BPR)*M_cruise) / sqrt((1+0.82*BPR)*G) ...
	+ X * p/p_sl * (0.23+0.19*sqrt(BPR)) * M_cruise^2);

% Installed thrust [N].
% Take a safety of around 20% more (according to T. Lambert).
% The plane actually need more thrust than for the cruise:
% we have to make turns, accelerations, maybe steep climbs, etc.
safety_thrust = 0.2;
F_installed = F_sl * (1 + safety_thrust);

% Thrust to weight ratio.
% The MTOW is conventionnaly expressed in [kg] and thus need to be
% converted in [N].
MTOW_force = MTOW * g;
TW_ratio = F_installed / MTOW_force;

% Uninstalled thrust [N].
% Lab geometrical conditions: assume optimistic configuration.
% We sould take into account the difference between the lab
% configuration and the actual placement of the engine on the plane:
% things like air duct and influence of wing/fuselage decrease the value
% of the installed thrust. Installed thrust is typically 4-8% lower than
% uninstalled thrust (slide 32, lesson 2).
install_loss = 8e-2;
F_uninstalled = F_installed / (1 - install_loss);

% Sanity check: verify that the engine is powerful enough.
assert(engine.Thrust >= F_uninstalled, ...
	"Thrust for the selected engine is too low");

%% Fuel weight estimation.
% Teams>General>Books>Raymer (page 149)

% Uninstalled SFC to cruise SFC.
% Teams>Propulsion>Files>Sources>Bartel-SFC_Calculations.pdf (page 12)
% Teams>Propulsion>Files>Sources>Nihad-Turbofans.pdf (page 33)
n = engine.n;
SFC_cruise = engine.SFC * sqrt(T/T_sl) * (1 + M_cruise)^n;
SFC_loiter = engine.SFC * sqrt(T/T_sl) * (1 + M_loiter)^n;

% Take into account takeoff and landing through a safety factor.
% According to T. Lambert, no need to worry too much about takeoff and
% landing at this design stage.
safety_takeoff_landing = 2e-2;
% Take into account the fact that the engine can have
% poorer-than-nominal fuel consumption.
% Teams>General>Books>Raymer (page 150)
safety_bad_engine = 5e-2;
% Take into account the fact that some fuel can be trapped in the fuel
% tank.
safety_trapped_fuel = 1e-2;

% Estimation of the thrust needed at loiter [N].
F_loiter = 0.5 * rho * V_loiter^2 * S * Cd;

% Summing amout of fuel required for all the mission steps [kg].
W_fuel = ...
	  (1+safety_takeoff_landing) ...
	* (1+safety_bad_engine+safety_trapped_fuel) ...
	* (  F_loiter * SFC_loiter * t_loiter1 ...
	   + F_loiter * SFC_loiter * t_loiter2 ...
	   + F_cruise * SFC_cruise * (eg_range/V_cruise) ...
	   + F_cruise * SFC_cruise * (ig_range/V_cruise));

% Fuel weight ratio.
FW_ratio = W_fuel / MTOW;

% Rough idea of the jet-A fuel density [kg/m³].
% Wikipedia>"Jet fuel"
rho_fuel = 820;

% Fuel tank volume [m³].
vol_tank = W_fuel / rho_fuel;
