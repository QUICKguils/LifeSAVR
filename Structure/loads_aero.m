function loads_aero()
% LOADS_AERO  Aerodynamic loads.
%
% This function computes the aerodynamic loads for all the critical
% points of the flight envelope.

%% Imports

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Constants and aircraft data.
C = load(fullfile(file_dir, "../constants.mat"));
D = load(fullfile(file_dir, "../data.mat"));

% Utilities.
addpath(genpath(fullfile(file_dir, "../Utils")));

%% Data.
theta_dd_add = 1.0472;  % Additional pitch acceleration [°/s²].
psi_max      = 15;      % Maximum yaw angle allowed [°].

end
