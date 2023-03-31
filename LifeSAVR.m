% TODO:
% - Finish to implement `opts`.

function LifeSAVR(opts)
% LIFESAVR  triggers all the code of the project.
%
% Parameter:
%	opts: char {'p', 'w'}, optional
%		'p' -> Enable plots creation.
%		'w' -> Write plotting data in external file.

%% Set path and global MAT files

% Find the root directory of the project.
root_dir = fileparts(mfilename('fullpath'));

% Add resursively sub-directories in the Matlab path.
addpath(genpath(fullfile(root_dir, "Utils")));
addpath(genpath(fullfile(root_dir, "Propulsion")));
addpath(genpath(fullfile(root_dir, "Structure")));

% Initialize MAT files.
constants();
data();

%% Options setting

% Option defaults: generate the plots, in ISoU.
if ~nargin
	opts = 'p';
end

%% Execute the code
% TODO: good idea to tell which function write in which MAT files.

% Load project constants.
C = load(fullfile(root_dir, "constants.mat"));

% Propulsion.
propulsion();
% pr_diagram();

% Structure.
flight_envelope(C.h_cr, opts);
% loads_aero();
% loads_structural();

% Debug print.
disp("LifeSAVR execution: OK");

end