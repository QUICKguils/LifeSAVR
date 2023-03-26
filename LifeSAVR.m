% TODO:
% - Add flexibility options.

function LifeSAVR(opts)
% LIFESAVR  triggers all the code of the project.
%
% Parameter:
%	opts: char {'p', 'w', 'i'}, optional
%		'p' -> Enable plots creation.
%		'w' -> Write data in external file.
%		'i' -> use the imperial system of units (abbr. ISoU).

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

% Determine if the imperial system of units is desired.
% TODO: see if it relly needs to be checked in this scope.
if contains(opts, 'i')
	ISoU = true;
else
	ISoU = false;
end

%% Execute the code
% TODO: good idea to tell which function write in  which MAT files.

% Propulsion.
propulsion();
% pr_diagram();
% 
% % Structure.
% flight_envelope();
% loads_aero();
% loads_structural();

% Debug print.
disp("exc OK");

end