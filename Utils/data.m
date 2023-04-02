function data
% DATA  Aircraft data used throughout the project.
%
% This function defines the aicrcraft data used throughout the project,
% and write them in `data.mat`.
%
% NOTE:
% This file should not exist. It should be removed, once all the matlab
% code is integrated and write itself in data.mat.

%% Imports

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Constants.
C = load(fullfile(file_dir, "../constants.mat"));

%% Wing

%Wing planform.
Wing.AR       = 8;         % Aspect ratio.
Wing.taper    = 0.35;      % Taper ratio.
Wing.surf     = 7.7339;    % Surface area [m²].
Wing.span     = 7.8658;    % Span [m].
Wing.c_root   = 1.4566;    % Chord at root [m].
Wing.c_tip    = 0.5098;    % Chord at tip [m].
Wing.smc      = 0.9832;    % Standard mean chord [m].
Wing.mac      = 1.0592;    % Mean Aerodynamic Chord [m].
Wing.y_ac     = 1.6508;    % Spanwise position of MAC [m].
Wing.twist    = -1;        % Aerodynamic twist (washout) [°].
Wing.sweep    = 24;        % Quarter-chord sweep angle [°].
Wing.mass     = 201.1955;  % Dry mass of wing [kg].
Wing.Vol_fuel = 0.39812;   % Volume of fuel in wing [m³].
Wing.cg_mac   = 30;        % Position of the CG at MAC w.r.t. the LE [%].
Wing.S_wet    = 15.9318;   % Wetted surface area [m²].
Wing.e        = 0.6164;    % Oswald's span efficiency.

% Wing Lift characteristics.
Wing.CL_alpha      = 5.6909;   % Slope of the lift coefficient [1/rad].
Wing.aoa_zerolift  = -1.7672;  % Zero-lift angle of attack at root [°].
Wing.aoa_midcruise = 0.52704;  % Angle of attack at root at mid-cruise [°].
Wing.aoi           = 0.9468;   % Angle of incidence [°].
Wing.CL_cr         = 0.29999;  % Refined Lift coefficient at cruise (MTOW) [-].

% Stall behaviour.
Wing.CL_max               = 1.4927;   % Max lift coefficient in clean configuration [-].
Wing.CL_max_takeoff       = 2.0164;   % Max lift coefficient at takeoff [-].
Wing.CL_max_landing       = 2.1458;   % Max lift coefficient at landing [-].
Wing.V_s_takeoff          = 57.7873;  % Stall velocity in takeoff configuration (MTOW) [m/s].
Wing.V_s_landing          = 40.3656;  % Stall velocity in landing configuration [m/s].
Wing.V_s_landing_mtow     = 56.0196;  % Stall velocity in landing configuration (MTOW) [m/s].
Wing.delta_cl_max_takeoff = 0.5237;   % Delta max lift coeff at takeoff.
Wing.delta_cl_max_landing = 0.6531;   % Delta max lift coeff at landing.

% Airfoil profile.
Wing.airfoil.name         = "NASA SC(2)-0412";
Wing.airfoil.tc_ratio     = 0.12;     % Thickness to chord ratio [-].
Wing.airfoil.cl_alpha     = 10.2106;  % Slope of the lift coefficient [1/rad].
Wing.airfoil.cl_max       = 1.72;     % Max lift coefficient (incompressible) [-].
Wing.airfoil.aoa_zerolift = -2.022;   % Zero-lift angle of attack of the airfoil [°].
Wing.airfoil.M_dd         = 0.786;    % Drag divergence mach number [-].
Wing.airfoil.cm           = -0.0778;  % Pitching moment coefficient (incompressible) [-].

% Flaps planform: plain flaps.
% Begins at the end of fuselage width + 2% of span
% Ends at 60% of span
% 25% of the airfoil chord
Wing.flaps.defl_takeoff = 20; % Flap deflection at takeoff [°].
Wing.flaps.defl_landing = 40; % Flap deflection at landing [°].

% Ailerons planform.
% Begins at 60% of span
% Ends at 90% of span
% 25% of the airfoil chord

%% Horizontal tail

HT.Vol_coeff        = 0.8;        % Volume coefficient.
HT.AR               = 4.25;       % Aspect ratio.
HT.taper            = 0.415;      % Taper ratio.
HT.sweep            = 29;         % Quarter-chord sweep angle [°].
HT.lever            = 3.58;       % Lever arm [m].
HT.surf             = 1.58;       % Surface [m²].
HT.visible_halfspan = 1.33;       % Individual visible span [m].
HT.c_root           = 0.84;       % Chord at root [m].
HT.c_tip            = 0.35;       % Chord at tip [m].
HT.mac              = 0.63;       % Mean aerodynamic chord [m].
HT.span             = 3.56;       % Total span [m].
HT.hidden           = 5;          % Percentage of hidden tail.
HT.W                = 15.6;       % Weight [kg].
HT.deda             = 0.45;       % Downwash slope [1/rad].
HT.eps              = 1.8;        % Downwash.
HT.CL               = -0.06;      % Lift coefficient.
HT.CL_alpha         = 1.8;        % Slope of the lift coefficient [1/rad].
HT.aoi              = -0.15;      % Angle of incidence [°].
HT.x_AC             = 0.19;       % X-position of the aerodynamic center [m].
HT.y                = 0.29;       % Distance from root to MAC [m].
HT.elevator.c       = 0.22;       % Elevator chord [m].
HT.elevator.L       = 1.2;        % Elevator span [m]
HT.elevator.sweep   = 23.8;       % Elevator sweep angle [°].
HT.height           = -0.2;       % Distance between wing and tail [m].
HT.CL_max           = -0.21;      % Max lift coefficient.
HT.airfoil.name     = "SC 0010";  % Airfoil name.

%% Vertical tail

VT.Vol_coeff     = 0.07;       % Volume coefficient.
VT.AR            = 1.5;        % Aspect ratio.
VT.taper         = 0.45;       % Taper ratio.
VT.sweep         = 24;         % Quarter-chord sweep angle [°].
VT.surf          = 1.2;        % Surface [m²].
VT.span          = 1.37;       % Span [m].
VT.c_root        = 1.2;        % Chord at root [m].
VT.c_tip         = 0.54;       % Chord at tip [m].
VT.mac           = 0.92;       % Mean aerodynamic chord [m].
VT.y             = 0.6;        % Distance from root to MAC [m].
VT.rudder.c      = 0.32;       % Rudder chord [m].
VT.rudder.length = 1.24;       % Rudder length [m].
VT.rudder.surf   = 0.4;        % Rudder surface [m²].
VT.rudder.sweep  = 13.7;       % Rudder sweep angle [°].
VT.W             = 37;         % Weight [kg].
VT.airfoil.name  = "SC 0010";  % Airfoil name.

%% Propulsion

% Propu structure, generated by Propulsion/propulsion.m.

%% Flight envelope

% FE structure, generated by Structure/flight_envelope.m

%% Components

Comp = table();
%                                                             COG [m]
%                 Name                    Mass [lb]      X       Y       Z
%    -----------------------------        ---------   -------|-------|------- 
Comp("Fuselage",                    :) = {1233.8,    [3.6,    0.0000, 0.0000]};
Comp("Main landing gear",           :) = {381.4,     [4.7,    0.0000, 0.0000]};
Comp("Nose landing gear",           :) = {80.2891,   [1.5,    0.0000, 0.0000]};
Comp("Engine mount",                :) = {6.7189,    [7.13,   0.0000, 0.0000]};
Comp("Engine section",              :) = {4.782,     [7.13,   0.0000, 0.0000]};
Comp("Air induction system",        :) = {333.1425,  [4.48,   0.0000, 0.0000]};
Comp("Tail pipe",                   :) = {4.7075,    [7.75,   0.0000, 0.0000]};
Comp("Starter",                     :) = {12.65,     [6.48,   0.0000, 0.0000]};
Comp("Fuel system and empty tanks", :) = {11.087,    [3.6028, 0.0000, 0.0000]};
Comp("Instruments",                 :) = {83.49,     [3.5,    0.0000, 0.0000]};
Comp("Hydraulics",                  :) = {171.7485,  [4.5,    0.0000, 0.0000]};
Comp("Furnishings",                 :) = {108.8,     [1.5,    0.0000, 0.0000]};
Comp("Handling gear",               :) = {0.6596,    [1.5,    0.0000, 0.0000]};
Comp("Vertical tail",               :) = {29.63,     [7.3595, 0.0000, 0.0000]};
Comp("Horizontal tail",             :) = {14.68,     [7.6576, 0.0000, 0.0000]};
Comp("Wings",                       :) = {414.83,    [4.1228, 0.0000, 0.0000]};
Comp("Fuel Wing",                   :) = {669.33908, [4.1228, 0.0000, 0.0000]};
Comp("Fuel fuselage",               :) = {3078.5149, [4.1228, 0.0000, 0.0000]};
Comp("Payload",                     :) = {300,       [4.1228, 0.0000, 0.0000]};
Comp("Engine",                      :) = {670.21,    [7.13,   0.0000, 0.0000]};
Comp("Flight control unit",         :) = {20,        [3,      0.0000, 0.0000]};
Comp("System control unit",         :) = {18,        [3,      0.0000, 0.0000]};
Comp("Flight transduces",           :) = {3,         [4.4202, 0.0000, 0.0000]};
Comp("Power engine control unit",   :) = {7,         [7.13,   0.0000, 0.0000]};
Comp("Power converter assembly",    :) = {4,         [7.13,   0.0000, 0.0000]};
Comp("GPS receiver (anti-jam)",     :) = {25,        [1.5,    0.0000, 0.0000]};
Comp("GPS antenna",                 :) = {6,         [1.5,    0.0000, 0.0000]};
Comp("Battery (non-propulsion)",    :) = {13,        [1.5,    0.0000, 0.0000]};
Comp("Cooling and pressurization",  :) = {5,         [2.6,    0.0000, 0.0000]};
Comp("Air data probe",              :) = {1,         [4.47,   0.0000, 0.0000]};
Comp("Anti-ice",                    :) = {2,         [1,      0.0000, 0.0000]};
Comp("Thermal management",          :) = {23,        [2.5,    0.0000, 0.0000]};
Comp("Miscelectrical wiring",       :) = {5,         [2,      0.0000, 0.0000]};
Comp("Basic satellite radio",       :) = {34,        [2.5,    0.0000, 0.0000]};
Comp("Data handling and linking",   :) = {24,        [2.5,    0.0000, 0.0000]};
Comp("Flight termination system",   :) = {10,        [2,      0.0000, 0.0000]};
Comp("Navigation fusion processor", :) = {9,         [2,      0.0000, 0.0000]};
Comp("Transponders",                :) = {15,        [1.5,    0.0000, 0.0000]};
Comp("Radar",                       :) = {63,        [0.75,   0.0000, 0.0000]};
Comp("Nose EO/IR/LIDAR",            :) = {100,       [1,      0.0000, 0.0000]};
Comp("Aft EO/IR/LIDAR",             :) = {50,        [6.5,    0.0000, 0.0000]};
Comp("Sensor growth",               :) = {15,        [1.5,    0.0000, 0.0000]};
Comp.Properties.VariableNames = {'Mass', 'COG'};

% Convert weights from lb to kg.
Comp.Mass = array2table(Comp.Mass .* C.lb2kg);

%% Plane

Plane.MTOW  = sum(Comp.Mass);                         % MTOW [kg].
Plane.COG = sum(Comp.Mass .* Comp.COG) / Plane.MTOW;  % COG [m].
Plane.CD_0  = 0.017;                                  % Zero lift drag coefficient.

%% Write in data.mat

% Save data in data.mat, which lies in the root directory.
save(fullfile(file_dir, "../data.mat"), "Wing", "HT", "VT", "Comp", "Plane");

end