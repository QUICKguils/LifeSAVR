function constants
% CONSTANTS  Constant quantities used thourghout the project.
%
% This function defines the constants used throughout the project, and
% write them in `constants.mat`.

clear;

%% Conversion factors.

in2m     = 2.54;         % Inch            ->  meter.
ft2m     = 0.3048;       % Foot            ->  meter.
mi2m     = 1609.344;     % Mile            ->  meter.
nmi2m    = 1852;         % Nautical mile   ->  meter.
mph2kmh  = 1609.344;     % Miles per hour  ->  kilometers per hour.
mph2ms   = 0.44704;      % Miles per hour  ->  meters per second.
kn2kmh   = 1.852;        % Knot            ->  kilometers per hour.
kn2ms    = 0.514444;     % Knot            ->  meters per second.
slug2kg  = 14.59390;     % Slug            ->  kilogram.
lb2kg    = 0.45359237;   % Pound mass      ->  kilogram.
lbf2N    = 4.448222;     % Pound force     ->  newton.
usgal2m3 = 3.785411784;  % US gallon       ->  cubic meter.


%% Physical constants.

g = 9.80665;                          % Gravitational acceleration [m/s²].
[rho_sl, p_sl, T_sl, a_sl] = ISA(0);  % Fluid properties at sea level. See Utils/ISA.m

%% Key performance parameters.
% Note that only the objective values are encoded.

% KPP_01 : irrelevant for our mission.
range_ingress      = 850 * nmi2m;  % KPP_02 Minimum ingress range [m].
M_ingress          = 0.86;         % KPP_03 Minimum ingress Mach number.
t_loiter_search    = 5 * 3600;     % KPP_04 Minimum loiter time to search for victims [s].
range_egress       = 2e3 * nmi2m;  % KPP_05 Minimum egress range [m].
M_egress           = 0.86;         % KPP_06 Minimum egress Mach number.
L_runway_landing   = 1500 * ft2m;  % KPP_07 Maximum landing rundway length [m].
M_dash             = 0.9;          % KPP_08 Minimum dash speed Mach number.
P_elec             = 24e3;         % KPP_09 Minimum available cruise electrical power [W].
t_loiter_landing   = 45 * 60;      % KPP_10 Minimum loiter time at landing site [s].
W_payload          = 300 * lb2kg;  % KPP_11 Minimum mission payload weight [kg].
n_ff               = 6;            % KPP_12 Minimum free flight loads.
t_swap             = 5 * 60;       % KPP_13 Maximum mission payload swap time [s].
t_repl_fluids      = 5 * 60;       % KPP_14 Maximum time to replenish operanting fluids [s].
t_repl_consumables = 5 * 60;       % KPP_15 Maximum time to replenish consumables [s].
N_lifespan         = 50;           % KPP_16 Minimum number of sorties before EOL of the plane.
L_runway_takeoff   = 1500 * ft2m;  % KPP_17 Maximum takeoff runway length [m].
n_up               = 3;            % KPP_18 Upward flight limit load factor.
n_down             = -1.5;         % KPP_19 Downward flight limit load factor.
vd_landing_lw      = 15 * ft2m;    % KPP_20 Landing maxiumum vertical descent at max. landing weight [m/s].
vd_landing_mtow    = 10 * ft2m;    % KPP_21 Landing maxiumum vertical descent at MTOW [m/s].
K_safety           = 1.5;          % KPP_22 Structural safety factor.

% Quantities directly derived from the KPPs.

h_cr     = 30e3 * ft2m;                  % Cruise altitude [m].
[rho_cr, p_cr, T_cr, a_cr] = ISA(h_cr);  % Fluid properties at cruise altitude. See Utils/ISA.m
M_cr     = M_ingress;                    % There's actually no difference bw ingress and egress.
V_cr     = M_cr * a_cr;                  % TAS at cruise [m/s].
M_loiter = 0.7;                          % Mach number during loiter.
V_loiter = M_loiter * a_cr;              % TAS at loiter [m/s].
V_dash   = M_dash * a_cr;                % TAS at dash   [m/s].

%% Save data into constants.mat

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Save data in constants.mat, which lies in the root directory.
save(fullfile(file_dir, "../constants.mat"));

end
