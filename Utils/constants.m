function constants
% CONSTANTS  Constant quantities used thourghout the project.
%
% This function defines the constants used throughout the project, and
% write them in `constants.mat`.

%% Conversion factors.

C.in2m     = 2.54;         % Inch            ->  meter.
C.ft2m     = 0.3048;       % Foot            ->  meter.
C.mi2m     = 1609.344;     % Mile            ->  meter.
C.nmi2m    = 1852;         % Nautical mile   ->  meter.
C.mph2kmh  = 1609.344;     % Miles per hour  ->  kilometers per hour.
C.mph2ms   = 0.44704;      % Miles per hour  ->  meters per second.
C.kn2kmh   = 1.852;        % Knot            ->  kilometers per hour.
C.kn2ms    = 0.514444;     % Knot            ->  meters per second.
C.slug2kg  = 14.59390;     % Slug            ->  kilogram.
C.lb2kg    = 0.45359237;   % Pound mass      ->  kilogram.
C.lbf2N    = 4.448222;     % Pound force     ->  newton.
C.usgal2m3 = 3.785411784;  % US gallon       ->  cubic meter.


%% Physical constants.

C.g = 9.80665;  % Gravitational acceleration [m/s²].
[C.rho_sl, C.p_sl, C.T_sl, C.a_sl] = ISA(0);  % Fluid properties at sea level. See Utils/ISA.m

%% Key performance parameters.
% Note that only the objective values are encoded.

% KPP_01 : irrelevant for our mission.
C.range_ingress      = 850 * C.nmi2m;  % KPP_02 Minimum ingress range [m].
C.M_ingress          = 0.86;           % KPP_03 Minimum ingress Mach number.
C.t_loiter_search    = 5 * 3600;       % KPP_04 Minimum loiter time to search for victims [s].
C.range_egress       = 2e3 * C.nmi2m;  % KPP_05 Minimum egress range [m].
C.M_egress           = 0.86;           % KPP_06 Minimum egress Mach number.
C.L_runway_landing   = 1500 * C.ft2m;  % KPP_07 Maximum landing rundway length [m].
C.M_dash             = 0.9;            % KPP_08 Minimum dash speed Mach number.
C.P_elec             = 24e3;           % KPP_09 Minimum available cruise electrical power [W].
C.t_loiter_landing   = 45 * 60;        % KPP_10 Minimum loiter time at landing site [s].
C.W_payload          = 300 * C.lb2kg;  % KPP_11 Minimum mission payload weight [kg].
C.n_ff               = 6;              % KPP_12 Minimum free flight loads.
C.t_swap             = 5 * 60;         % KPP_13 Maximum mission payload swap time [s].
C.t_repl_fluids      = 5 * 60;         % KPP_14 Maximum time to replenish operanting fluids [s].
C.t_repl_consumables = 5 * 60;         % KPP_15 Maximum time to replenish consumables [s].
C.N_lifespan         = 50;             % KPP_16 Minimum number of sorties before EOL of the plane.
C.L_runway_takeoff   = 1500 * C.ft2m;  % KPP_17 Maximum takeoff runway length [m].
C.n_up               = 3;              % KPP_18 Upward flight limit load factor.
C.n_down             = -1.5;           % KPP_19 Downward flight limit load factor.
C.vd_landing_lw      = 15 * C.ft2m;    % KPP_20 Landing maxiumum vertical descent at max. landing weight [m/s].
C.vd_landing_mtow    = 10 * C.ft2m;    % KPP_21 Landing maxiumum vertical descent at MTOW [m/s].
C.K_safety           = 1.5;            % KPP_22 Structural safety factor.

% Quantities directly derived from the KPPs.

C.h_cr     = 30e3 * C.ft2m;        % Cruise altitude [m].
[C.rho_cr, C.p_cr, C.T_cr, C.a_cr] = ISA(C.h_cr);  % Fluid properties at cruise altitude. See Utils/ISA.m
C.M_cr     = C.M_ingress;          % There's actually no difference bw ingress and egress.
C.V_cr     = C.M_cr * C.a_cr;      % TAS at cruise [m/s].
C.M_loiter = 0.7;                  % Mach number during loiter.
C.V_loiter = C.M_loiter * C.a_cr;  % TAS at loiter [m/s].
C.V_dash   = C.M_dash * C.a_cr;    % TAS at dash   [m/s].

%% Save data into constants.mat

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Save data in constants.mat, which lies in the root directory.
save(fullfile(file_dir, "../constants.mat"), '-struct', "C");

end
