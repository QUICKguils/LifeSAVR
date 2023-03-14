function [rho, p, T, a] = ISA(h)
% ISA  International standard atmosphere.
%
% Returns useful flow properties for the given altitude, in meter.
%
% Parameter:
%	h: double
%	  Altitude, in meter.
% Return:
%	[rho, p, T, a]: 1x4 double
%	  Density, pressure, temperature and sound speed, in SI units.

% Constants and quantities at sea level.
T_sl   = 288.15;   % Temperature [K].
p_sl   = 101325;   % Pressure [Pa].
rho_sl = 1.225;    % Density of air [kg/m³].
g      = 9.80665;  % Gravitational acceleration [m/s²].
R      = 287.04;   % Molar gas constant of air [m²/(s²K)].
aLR    = -0.0065;  %
gamma  = 1.4;      % Heat capacity ratio of air.

% Temperature, pressure and density.
if h < 11e3
	T   = T_sl + aLR*h;
	p   = p_sl   * (T/T_sl)^(  -g/(aLR*R));
	rho = rho_sl * (T/T_sl)^(-1-g/(aLR*R));
else
	theta = 0.7519;
	if h <= 36089
		delta = 0.223 * exp(-0.0001578*(h-11000));
		sigma = 0.297 * exp(-0.0001578*(h-11000));
	else
		delta = 0.223 * exp(-0.00004811*(h-36089));
		sigma = 0.297 * exp(-0.00004811*(h-36089));
	end
	T   = T_sl * theta;
	p   = p_sl * delta;
	rho = rho_sl * sigma;
end

% Sound speed.
a = sqrt(gamma * R * T);
end