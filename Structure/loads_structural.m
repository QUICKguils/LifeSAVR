function loads_structural()
% LOADS_STRUCT  Structural loads.

% This function computes the structural loads for all the critical
% points of the flight envelope, along all the wing and fuselage cross
% sections.

%% Imports

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Constants and aircraft data.
C = load(fullfile(file_dir, "../constants.mat"));
D = load(fullfile(file_dir, "../data.mat"));

% Utilities.
addpath(genpath(fullfile(file_dir, "../Utils")));

end
