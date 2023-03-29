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

% Air properties at the desired altitude.
[rho, ~, T, a] = ISA(h);

% From the statement.
theta_dd_add = 1.0472;  % Additional pitch acceleration [rad/s²].
psi_max      = 15;      % Maximum yaw angle allowed [°].

% Critical points of the flight envelope.
V = Envelope.CP.n;
n = Envelope.CP.EAS;

%% Aircraft loads.

result_alpha = zeros(1,length(V));
result_L_w   = zeros(1,length(V));
result_L_t   = zeros(1,length(V));

for i = 1:numel(V)

    % Thrust [N].
	% WARN: check if TAS or EAS.
	% Conversion of true airspeed (TAS) to equivalent airspeed (EAS).
	TAS2EAS = sqrt(rho/C.rho_sl);
	M = TAS * TAS2EAS/a;
    T(i) = thrust_sl2alt(V(i), h);
    % Weight [N].
    nW(i) = D.Plane.MTOW * C.g * n(i);
    % Drag [N].
	% TODO: must depend on the envelope point.
	CD_wing = Wing.CD_0 + Wing.CL_cr^2 / (Wing.e * pi * Wing.AR);
	% C_D_wing = max([c_D_cruise_loiter, c_D_cruise_in,  c_D_cruise_eg]);
    D_w(i) = 0.5 * CD_wing * V(i)^2 * rho * D.Wing.S;
	% TODO: modify when we have more precise values.
    CD_body = 0.015;
    D_b(i)  = 0.5 * CD_body * V(i)^2 * rho * 2*pi*width_f^2/4;
    % Pitching moment.
    C_M_0 = -0.1;
    K_n   = 0.09; % TODO: Check if this value is OK.
	% This formula does not apply to our geometry.
    % C_M = C_M_0 - K_n * D.Wing.CL_max - 0.0015 * psi_max;
	C_M = D.Wing.airfoil.C_m;
    M(i) = 0.5 * C_M * V(i)^2 * rho * D.Wing.S * D.Wing.smc;

    I_theta = 3000; % A calculer avec le CAD
    c_l_alpha_plane = 6.38;
    alpha_L0_f = -0.0539;

	% TODO: change once we have the z_pos of the required components
	% (wait for the CAD).
	% NOTE: L_T c'est le P des formules.
    syms alpha_var L_w L_t
    eqI   = c_l_alpha_plane * (alpha_var - alpha_L0_f) - L_w / (0.5 * C.rho_cr * D.Wing.S * (V(i)^2)); %Lift of the wing
    eqII  = abs(c_g - pos(1)) * L_w - abs(z_cg - z_pos(8)) * T(i) + z_cg * D_b(i) - z_pos(1) * D_w(i) - abs(pos(4) - c_g) * L_t + M(i) - I_theta * theta_dd_add; % Moments about COG
    eqIII = L_w + L_t - nW(i) + T(i) * sin(alpha_var-alpha_L0_f); % Vertical equilibrium

    [alpha_var, L_w, L_t]=vpasolve(eqI, eqII, eqIII, alpha_var, L_w, L_t);

    result_alpha(i) = double(rad2deg(alpha_var)); %[deg]
    result_L_w(i)   = double(L_w);
    result_L_t(i)   = double(L_t);

	% Lift curve slope of the vertical tail.
	a_VT = 5.5 * AR_VT / (AR_VT + 2);
	% Lift of the vertical tail.
	F_fin(i) = 0.5 * rho * V(i)^2 * a_VT* S_VT * psi_max;

	% NOTE: no need to take into account the M_tail.

	% Fuselage bending moment.
    M_fus(i) = F_fin(i) * (z_cg - z_pos(3));
end

%% Fuselage loads.

x_rear_fus = pos(1) + 3/4*c_root_w;
l_rear_fus = L_f-x_rear_fus;
W_rear_fus = l_rear_fus/L_f * W(5)*g_SI;

l_tail_fin = pos(3)-x_rear_fus;
W_tail_fin = (W(3) + W(4))*g_SI;

l_engine = pos(8) - x_rear_fus;
W_engine = W(8)*g_SI;

x_dd = L_f - x_rear_fus;
x_bb = L_f - (pos(8) - L_eng/2); %beginning of the engine
x_cc = x_bb/2;

D = 0.8; %D is constant over the entire rear fuselage length

A_skin = pi*D*(L_f-x_rear_fus);

weight_m = W_rear_fus* D* pi/A_skin; %the same for each section because the diameter is constant

for i = 1:length(V)
    S_FA(i) = n_points(i)*(W_rear_fus + W_tail_fin + W_engine);
    B_MA(i) = n_points(i)*cosd(result_alpha(i)-0.08)*(l_rear_fus/2*W_rear_fus + l_tail_fin*W_tail_fin * l_engine*W_engine);

    T_yA(i) = - F_fin(i);
    T_zA(i) = - (S_FA(i) - result_L_t(i)) * cosd(result_alpha(i)-0.08);
    M_yA(i) =   B_MA(i) - result_L_t(i) * l_tail_fin * cosd(result_alpha(i)-0.08);
    M_zA(i) =   F_fin(i) * l_tail_fin;
    M_xA(i) = - M_tail(i) - F_fin(i) * abs(z_pos(3));

    S_FC(i) = n_points(i) * ((l_rear_fus-x_cc) / L_f * W_rear_fus + W_tail_fin + W_engine);
    B_MC(i) = n_points(i)*cosd(result_alpha(i)-0.08)*((l_rear_fus-x_cc)^2/(2*L_f)*W_rear_fus + (l_tail_fin-x_cc)*W_tail_fin * (l_engine-x_cc)*W_engine);

    T_yC(i) = - F_fin(i);
    T_zC(i) = - (S_FC(i) - result_L_t(i))*cosd(result_alpha(i)-0.08);
    M_yC(i) =   B_MC(i) - result_L_t(i)*(l_tail_fin-x_cc)*cosd(result_alpha(i)-0.08);
    M_zC(i) =   F_fin(i) * (l_tail_fin - x_cc);
    M_xC(i) = - M_tail(i) - F_fin(i) * abs(z_pos(3));

    S_FB(i) = n_points(i)*((l_rear_fus-x_bb)/L_f*W_rear_fus + W_tail_fin + W_engine);
    B_MB(i) = n_points(i)*cosd(result_alpha(i)-0.08)*((l_rear_fus-x_bb)^2/(2*L_f)*W_rear_fus + (l_tail_fin-x_bb)*W_tail_fin * (l_engine-x_bb)*W_engine);

    T_yB(i) = - F_fin(i);
    T_zB(i) = - (S_FB(i) - result_L_t(i))*cosd(result_alpha(i)-0.08);
    M_yB(i) =   B_MC(i) - result_L_t(i)*(l_tail_fin-x_bb)*cosd(result_alpha(i)-0.08);
    M_zB(i) =   F_fin(i) * (l_tail_fin-x_bb);
    M_xB(i) = - M_tail(i) - F_fin(i)*abs(z_pos(3));
end

end
