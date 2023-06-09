% PROCEDURE:
% - Get an estimation of the thrust needed for the cruise speed.
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
%   maybe add something like 2%).
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
% - Check the 20% and the 8% factor: 20% seems too high. Provide refs.
% - Think more about the air duct. What size, how to beneficiate from
%   the BL of the nose/fuselage. See lesson 5 (Koen).
% - Carry out a small tradeoff study for the engine: BRP and SFC vs
%   engine diameter (fuselage drag).

function propulsion
% PROPULSION  Calculations related to the propulsion.
%
% Save:
%   Propu: struct
%     Holds all the relevant data generated by this function.

%% Imports

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Constants and aircraft data.
C = load(fullfile(file_dir, "../constants.mat"));
D = load(fullfile(file_dir, "../data.mat"));

% Utilities.
addpath(genpath(fullfile(file_dir, "../Utils")));

%% Engine selection

% Selected engine data.
% Thrust, SFC and BPR data are gathered by looking at either the
% manufacturer's data sheets or the engine certificate (FAA or EASA).
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
engine = engine_table("FJ44-4A", :);

%% Estimation of the thrust needed for cruise speed
% Slides 27 and following, lesson 2.

% Estimation of e is now calculated in the wing part.
% % Estimation of the Oswald's efficiency factor.
% % Teams>General>Books>Raymer (page 444)
% e_straight = 1.78 * (1 - 0.045*AR^0.68)                     - 0.64;
% e_swept    = 4.61 * (1 - 0.045*AR^0.68) * cosd(Lambda)^0.15 - 3.1;
% % Linear interpolation for a sweep angle of Lambda.
% e = e_straight + (e_swept - e_straight)/30 * D.Wing.sweep;

% Estimation of Cd in now calculated in the drag part.
% % Estimation of the drag coefficient.
% Cd = D.Plane.CD_0 + D.Wing.CL_cr^2 / (D.Wing.e * pi * D.Wing.AR);

% Thrust at cruise must equal the drag [N].
F_cr = 0.5 * C.rho_cr * C.V_cr^2 * D.Wing.surf * D.Plane.CD_cr;

%% Cruise thrust to uninstalled thrust
% The uninstalled thrust corresponds to the thrust measured in lab
% conditions by the manufacturers.
% WARN: conditions assumed below can depend on the manufacturer!

% Thrust at sea level, static (SLS) [N].
% Actually there is a causal dilemma: we want to estimate the equivalent
% uninstalled thrust to select a matching engine. But at the same time,
% we need to know which engine is used in order to evaluate this
% equivalent uninstalled thrust. More precisely, the equivalent thrust
% depend on the engine BPR.
% -> We can have a first guesstimate of 0.8 for the value of G. This
% corresponds to a quite small BPR engine.
% Teams>Propulsion>Files>Sources>Bartel-SFC_Calculations.pdf (page 8)
SLS2cr = thrust_SLSconv(C.V_cr, C.h_cr, engine.BPR, engine.G);
F_sls = F_cr / SLS2cr;

% Installed thrust [N].
% Take a safety of around 15% more (according to T. Lambert).
% The plane actually need more thrust than for the cruise:
% we have to make turns, accelerations, maybe steep climbs, etc.
safety_thrust = 0.15;
F_installed = F_sls * (1 + safety_thrust);

% Thrust to weight ratio.
% The MTOW is conventionally expressed in [kg] and thus need to be
% converted in [N].
TW_ratio = F_installed / (D.Plane.MTOW * C.g);

% Uninstalled thrust [N].
% Lab geometrical conditions: assume optimistic configuration.
% We should take into account the difference between the lab
% configuration and the actual placement of the engine on the plane:
% things like air duct and influence of wing/fuselage decrease the value
% of the installed thrust. Installed thrust is typically 4-8% lower than
% uninstalled thrust (slide 32, lesson 2).
install_loss = 8e-2;
F_uninstalled = F_installed / (1 - install_loss);

% Sanity check: verify that the engine is powerful enough.
assert(engine.Thrust >= F_uninstalled, ...
	"Thrust for the selected engine is too low");

%% Design cruise Mach and design dive Mach
% The design cruise Mach is defined as the Mach number obtained at
% maximum engine thrust, in cruising conditions.
% The design dive Mach can be calculated as 1.25 times the design Mach.
% Aircraft Structures>lesson 1>slides 7 and 11.

% Maximum thrust at cruise [N].
% We just have to revert steps of the preceding section, given the
% uninstalled thrust of the selected engine.
F_c = engine.Thrust * (1 - install_loss) * SLS2cr;

% Design cruise speed [m/s].
% Obtained through the thrust-drag equality.
V_c = sqrt(F_c / (0.5 * C.rho_cr * D.Wing.surf * D.Plane.CD_cr));

% Design cruise and dive mach numbers.
M_c = V_c / C.a_cr;
M_d = 1.25 * M_c;

%% Fuel weight estimation
% Teams>General>Books>Raymer (page 149)

% Uninstalled SFC to cruise SFC.
% Teams>Propulsion>Files>Sources>Bartel-SFC_Calculations.pdf (page 12)
% Teams>Propulsion>Files>Sources>Nihad-Turbofans.pdf (page 33)
n = engine.n;
SFC_cruise = engine.SFC * sqrt(C.T_cr/C.T_sl) * (1 + C.M_cr)^n;
SFC_loiter = engine.SFC * sqrt(C.T_cr/C.T_sl) * (1 + C.M_loiter)^n;

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
F_loiter = 0.5 * C.rho_cr * C.V_loiter^2 * D.Wing.surf * D.Plane.CD_loiter;

% Summing amout of fuel required for all the mission steps [kg].
W_fuel = ...
	  (1+safety_takeoff_landing) ...
	* (1+safety_bad_engine+safety_trapped_fuel) ...
	* (  F_loiter * SFC_loiter * C.t_loiter_search ...
	   + F_loiter * SFC_loiter * C.t_loiter_landing ...
	   + F_cr     * SFC_cruise * (C.range_egress/C.V_cr) ...
	   + F_cr     * SFC_cruise * (C.range_ingress/C.V_cr));

% Fuel weight ratio.
FW_ratio = W_fuel / D.Plane.MTOW;

% Rough idea of the jet-A fuel density [kg/m³].
% Wikipedia>"Jet fuel"
rho_fuel = 820;

% Fuel tank volume [m³].
vol_tank = W_fuel / rho_fuel;

%% Update values in data.mat

% Collect the relevant quantities to save, in a structure named `Propu`.
Propu.engine_table  = engine_table;
Propu.engine        = engine;
Propu.T_cr          = F_cr;
Propu.T_sls         = F_sls;
Propu.T_installed   = F_installed;
Propu.T_uninstalled = F_uninstalled;
Propu.TW_ratio      = TW_ratio;
Propu.T_c           = F_c;
Propu.M_c           = M_c;
Propu.M_d           = M_d;
Propu.SFC_cruise    = SFC_cruise;
Propu.SFC_loiter    = SFC_loiter;
Propu.W_fuel        = W_fuel;
Propu.FW_ratio      = FW_ratio;
Propu.vol_tank      = vol_tank;

% Save Propu in data.mat.
save(fullfile(file_dir, "../data.mat"), "Propu", "-append");

end