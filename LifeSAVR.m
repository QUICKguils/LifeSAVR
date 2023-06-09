function LifeSAVR(opts)
% LIFESAVR  triggers all the code of the project.
%
% Argument:
%	opts: char {'p', 'w', 'i'}, optional. Default is 'pi'.
%	  'p' -> Enable plots creation.
%	  'w' -> Write plotting data in external file.
%     'i' -> Plot and write data with imperial system of units (abbr. ISoU).

close all;

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
	opts = 'pi';
end

%% Execute the code
% TODO: good idea to tell which function write in which MAT files.

% Load project constants.
C = load(fullfile(root_dir, "constants.mat"));

% Propulsion.
propulsion();
pr_diagram('s', opts);

% Structure.
flight_envelope(C.h_cr, opts);
aero_loads();
struct_loads();
struct_stress(opts);

end