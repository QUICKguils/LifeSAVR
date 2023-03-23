% TODO:
% this file should be removed, once all the matlab code is integrated
% and write themselves in data.mat.
%
% Last update from teams: 23/03

clear;

%% Wing

%Wing planform.
Wing.AR     = 8;         % Aspect ratio.
Wing.lambda = 0.35;      % Taper ratio.
Wing.S      = 7.3686;    % Surface area [m²].
Wing.b      = 7.6778;    % Span [m].
Wing.c_r    = 1.4218;    % Chord at root [m].
Wing.c_t    = 0.4976;    % Chord at tip [m].
Wing.smc    = 0.9597;    % Standard mean chord [m].
Wing.mac    = 1.0338;    % Mean Aerodynamic Chord [m].
Wing.y_ac   = 1.6114;    % Spanwise position of MAC [m].
Wing.twist  = -1;        % Aerodynamic twist (washout) [°].
Wing.sweep  = 24;        % Sweep angle (1/4 chord) [°].
Wing.mass   = 188.96;    % Dry mass of wing [kg].
Wing.V_fuel = 0.370253;  % Volume of fuel in wing [m³].
Wing.cg_mac = 30;        % Position of the CG of the wing at MAC w.r.t. the LE [%].
Wing.S_wet  = 15.1793;   % Wetted surface area [m²].
Wing.e      = 0.6164;    % Oswald's span efficiency.

% Wing Lift characteristics.
Wing.CL_alpha      = 5.6909;   % Slope of the lift coefficient of the wing [1/rad].
Wing.aoa_zerolift  = -1.7672;  % Zero-lift angle of attack at root [°].
Wing.aoa_midcruise = 0.5044;   % Angle of attack at root at mid-cruise [°].
Wing.i             = 0.9241;   % Angle of incidence of the wing [°].
Wing.CL_cr         = 0.276;    % Refined Lift coefficient at cruise (MTOW) [-].

% Wing drag characteristics.
Wing.CD_0 = 0.017;  % Zero lift drag coefficient.

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

% Airfoil profile.
Wing.airfoil.name         = "NASA SC(2)-0412";
Wing.airfoil.tc_ratio     = 0.12;     % Thickness to chord ratio [-].
Wing.airfoil.cl_alpha     = 10.2106;  % Slope of the lift coefficient [1/rad].
Wing.airfoil.cl_max       = 1.72;     % Max lift coefficient (incompressible) [-].
Wing.airfoil.aoa_zerolift = -2.022;   % Zero-lift angle of attack of the airfoil [°].
Wing.airfoil.M_DD         = 0.786;    % Drag divergence mach number [-].
Wing.airfoil.C_m          = -0.0778;  % Pitching moment coefficient (incompressible) [-].

% Stall behaviour.
Wing.CL_max_clean         = 1.4927;   % Wing max lift coefficient in clean configuration [-].
Wing.CL_max_takeoff       = 2.0397;   % Wing max lift coefficient at takeoff [-].
Wing.CL_max_landing       = 2.1723;   % Wing max lift coefficient at landing [-].
Wing.Vs_takeoff           = 55.1174;  % Stall velocity in takeoff configuration (MTOW) [m/s].
Wing.Vs_landing           = 36.2724;  % Stall velocity in landing configuration [m/s].
Wing.Vs_landing_mtow      = 53.409;   % Stall velocity in landing configuration (MTOW) [m/s].
Wing.delta_cl_max_takeoff = 0.55;     % Delta max lift coeff at takeoff.
Wing.delta_cl_max_landing = 0.68;     % Delta max lift coeff at landing.

%% Tail



%% Write in data.mat

save("../data.mat", "-append")