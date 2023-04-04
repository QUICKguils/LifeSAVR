% Data of the material used

% TODO:
% - update the material propreties value when aurore finish the
%   material selection part

sigma_y   = 280e6;            % [Pa]
tau_y     = sigma_y/sqrt(3);  % [Pa]
s         = 1.5;              % safety factor (KPP.22)
sigma_max = sigma_y/s;        % [Pa]
tau_max   = tau_y/s;          % [Pa]

%% Fuselage design

n_stringers = 24;     % arbirary
l = 2*pi/n_stringers; % space between stringers

a = [1.3, 1.3, 1.3, 1.3];  %minor axis
b = [1, 1, 1, 1];           %major axis
a_aa = a(1);
a_bb = a(2);
a_cc = a(3);
a_dd = a(4);
b_aa = b(1);
b_bb = b(2);
b_cc = b(3);
b_dd = b(4);

% Coordinates of the stringers

for j = 1:length(a)-1
	for i = 0:n_stringers
		r(j) = sqrt(a(j)^2 * b(j)^2/( b(j)^2 * cos(i*l)^2 + a(j)^2 * sin(i*l)^2));  %TODO: not sure about the formula used.
		z_i(i+1,j) = r(j)*sin(-l*r(j)*(i));
		y_i(i+1,j) = r(j)*cos(-l*r(j)*(i));
	end
end

% TODO:
% - do we need to take the maximum value of My to find the B_min or it will be
%   found automatically at the end?

max_My_aa = find(max(abs(My_aa))); %the position at which we have the maximum moment on the vector we take the max value in order tu find the minimum value of B
My_aa_max = My_aa(max_My_aa);


max_My_bb = find(max(abs(My_bb)));
My_bb_max = My_bb(max_My_bb);

max_My_cc = find(max(abs(My_cc)));
My_cc_max = My_cc(max_My_cc);

max_Mz_aa = find(max(abs(Mz_aa)));
Mz_aa_max = Mz_aa(max_Mz_aa);

max_Mz_bb = find(max(abs(Mz_bb)));
Mz_bb_max = Mz_bb(max_Mz_bb);

max_Mz_cc = find(max(abs(Mz_cc)));
Mz_cc_max = Mz_cc(max_Mz_cc);

%-----inertia per unit boom------

I2B_yy_aa = sum(z_i(:,1).^2);  %Iy/B = sum(z)
I2B_zz_aa = sum(y_i(:,1).^2);  %Iz/B = sum(y)

I2B_yy_bb = sum(z_i(:,3).^2);
I2B_zz_bb = sum(y_i(:,3).^2);

I2B_yy_cc = sum(z_i(:,2).^2);
I2B_zz_cc = sum(y_i(:,2).^2);

%I_yz = 0; because we have symmetry

% Direct stress

B_sigma_aa = (My_aa_max*z_i(:,1)/I2B_yy_aa) - (Mz_aa_max*y_i(:,1)/I2B_zz_aa);
B_sigma_bb = (My_bb_max*z_i(:,3)/I2B_yy_bb) - (Mz_bb_max*y_i(:,3)/I2B_zz_bb);
B_sigma_cc = (My_cc_max*z_i(:,2)/I2B_yy_cc) - (Mz_cc_max*y_i(:,2)/I2B_zz_cc);

% minimum area of the stringers 

B_aa_min = max(abs(B_sigma_aa))/sigma_max;
B_bb_min = max(abs(B_sigma_bb))/sigma_max;
B_cc_min = max(abs(B_sigma_cc))/sigma_max;

% skin thickness symmetry is assumed about y and z axis if we neglect the air
% induction system

% TODO:
% - do we need to take the maximum value of My to find the B_min or it
%   will be found automatically at the end ?

max_Ty_aa = find(max(abs(Ty_aa)));
Ty_aa_max = Ty_aa(max_Ty_aa);

max_Ty_bb = find(max(abs(Ty_bb)));
Ty_bb_max = Ty_bb(max_Ty_bb);

max_Ty_cc = find(max(abs(Ty_cc)));
Ty_cc_max = Ty_cc(max_Ty_cc);

max_Tz_aa = find(max(abs(Tz_aa)));
Tz_aa_max = Tz_aa(max_Tz_aa);

max_Tz_bb = find(max(abs(Tz_bb)));
Tz_bb_max = Tz_bb(max_Tz_bb);

max_Tz_cc = find(max(abs(Tz_cc)));
Tz_cc_max = Tz_cc(max_Tz_cc);

max_Mx_aa = find(max(abs(Mx_aa)));
Mx_aa_max = Mx_aa(max_Mx_aa);

max_Mx_bb = find(max(abs(Mx_bb)));
Mx_bb_max = Mx_bb(max_Mx_bb);

max_Mx_cc = find(max(abs(Mx_cc)));
Mx_cc_max = Mx_cc(max_Mx_cc);

% shear due to T_z
q_z_aa(1) = -1/2*Tz_aa_max/I2B_yy_aa*z_i(1,1);
q_z_bb(1) = -1/2*Tz_bb_max/I2B_yy_bb*z_i(1,3);
q_z_cc(1) = -1/2*Tz_cc_max/I2B_yy_cc*z_i(1,2);

for i = 2:n_stringers
	q_z_aa(i) = q_z_aa(i-1) - Tz_aa_max/I2B_yy_aa*z_i(i,1);
	q_z_bb(i) = q_z_bb(i-1) - Tz_bb_max/I2B_yy_bb*z_i(i,3);
	q_z_cc(i) = q_z_cc(i-1) - Tz_cc_max/I2B_yy_cc*z_i(i,2);
end

% shear due to T_y
q_y_aa(1) = -1/2*Ty_aa_max/I2B_zz_aa*y_i(1,1);
q_y_bb(1) = -1/2*Ty_bb_max/I2B_zz_bb*y_i(1,3);
q_y_cc(1) = -1/2*Ty_cc_max/I2B_zz_cc*y_i(1,2);
for i = 2:n_stringers
	q_y_aa(i) = q_y_aa(i-1) - Ty_aa_max/I2B_zz_aa*y_i(i,1);
	q_y_bb(i) = q_y_bb(i-1) - Ty_bb_max/I2B_zz_bb*y_i(i,3);
	q_y_cc(i) = q_y_cc(i-1) - Ty_cc_max/I2B_zz_cc*y_i(i,2);

end

% TODO: wait for the CAD
A_aa = pi * D_aa^2 / 4;
A_bb = pi * D_bb^2 / 4;
A_cc = pi * D_cc^2 / 4;

% shear due to torque M_x
q_T_aa = Mx_aa_max/(2*A_aa); %A_aa = pi*D^2/4: is the area of section AA (given by the CAD)
q_T_bb = Mx_bb_max/(2*A_bb);
q_T_cc = Mx_cc_max/(2*A_cc);

% Total shear flow
q_tot_aa = q_T_aa + q_z_aa + q_y_aa;
q_tot_bb = q_T_bb + q_z_bb + q_y_bb;
q_tot_cc = q_T_cc + q_z_cc + q_y_cc;

% maximum shear flow
q_max_aa = max(abs(q_tot_aa));
q_max_bb = max(abs(q_tot_bb));
q_max_cc = max(abs(q_tot_cc));

% thickness
t_min_aa = q_max_aa/tau_max;
t_min_bb = q_max_bb/tau_max;
t_min_cc = q_max_cc/tau_max;

t_min = max([t_min_aa t_min_bb t_min_cc]);

% Rivets size
% not very interesting information (not mentionned in other reports)

% R_aa = (-Ty_aa/I2B_zz_aa*y_i(:,1))-(Tz_aa/I2B_yy_aa*z_i(:,1));
% R_bb = (-Ty_bb/I2B_zz_bb*y_i(:,3))-(Tz_bb/I2B_yy_bb*z_i(:,3));
% R_cc = (-Ty_cc/I2B_zz_cc*y_i(:,2))-(Tz_cc/I2B_yy_cc*z_i(:,2));

% R_max = max(abs([R_aa R_bb R_cc]));
