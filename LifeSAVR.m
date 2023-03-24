% TODO:
% - Add flexibility options.

function LifeSAVR()
% LIFESAVR  triggers all the code of the project.
%
% Parameter:
%	opts: char {'p', 'w', 'i'}, optional
%		'p' -> Enable plots creation.
%		'w' -> Write data in external file.
%		'i' -> use the imperial system of units (abbr. ISoU).

%% Update constants and data.

run Utils\constants.m
run Utils\data.m

%% Execute the code.

run Propulsion\propulsion.m  % -> data.m>Propu
run Propulsion\pr_diagram.m

run Structure\flight_envelope.m   % -> data.m>Struct
run Structure\loads_aero.m        % -> data.m>Struct
run Structure\loads_structural.m  % -> data.m>Struct

end