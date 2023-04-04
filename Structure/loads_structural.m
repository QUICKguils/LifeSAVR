%% Data
L_fuselage = 8;
L_engine = 1.74;

CP = D.FE.CP;

x_rear_fuselage = D.Comp{ "Wings", "COG"}(1) + 3/4 * D.Wing.c_root;
l_rear_fuselage = L_fuselage - x_rear_fuselage;
W_rear_fuselage = l_rear_fuselage * D.Comp{"Fuselage" , "Mass"} * C.g;

l_tail_fin = D.Comp{"Horizontal tail" , "COG"}(1) - x_rear_fuselage;
W_tail_fin = (D.Comp{"Vertical tail" , "Mass"} + D.Comp{"Horizontal tail", "Mass"}) * C.g;

l_engine = D.Comp{"Engine", "COG"}(1) - x_rear_fuselage;
W_engine = D.Comp{"Engine" , "Mass"} * C.g;

x_dd = L_fuselage - x_rear_fuselage;
x_bb = (D.Comp{"Engine", "COG"}(1) - L_engine/2) - x_rear_fuselage;  % beginning of the engine
x_cc = x_bb/2;

% D = [D_aa D_bb D_cc D_dd]; %[D_A D_B D_C D_D], D is constant over the entire rear fuselage length 
Di = [4, 5, 6, 7];
D_aa = Di(1);
D_bb = Di(2);
D_cc = Di(3);
D_dd = Di(4);

A_skin = pi*(max(Di)+min(Di))/2*l_rear_fuselage;

 W_maa = W_rear_fuselage * D_aa*pi/A_skin;
 W_mbb = W_rear_fuselage * D_bb*pi/A_skin;
 W_mcc = W_rear_fuselage * D_cc*pi/A_skin;
 W_mdd = W_rear_fuselage * D_dd*pi/A_skin;
 
 %% MNT

CP_num = height(CP);

SF_A = zeros(1, CP_num);
SF_B = zeros(1, CP_num);
SF_C = zeros(1, CP_num);

BM_A = zeros(1, CP_num);
BM_B = zeros(1, CP_num);
BM_C = zeros(1, CP_num);

Ty_aa = zeros(1, CP_num);
Tz_aa = zeros(1, CP_num);
My_aa = zeros(1, CP_num);
Mz_aa = zeros(1, CP_num);
Mx_aa = zeros(1, CP_num);

Ty_cc = zeros(1, CP_num);
Tz_cc = zeros(1, CP_num);
My_cc = zeros(1, CP_num);
Mz_cc = zeros(1, CP_num);
Mx_cc = zeros(1, CP_num);

Ty_bb = zeros(1, CP_num);
Tz_bb = zeros(1, CP_num);
My_bb = zeros(1, CP_num);
Mz_bb = zeros(1, CP_num);
Mx_bb = zeros(1, CP_num);

	result_alpha = D.AeroLoads.aoa;
	
	result_L_t = D.AeroLoads.P;
 for i = 1:CP_num
	n   = CP.n(i);
	F_fin = D.AeroLoads.F_fin(i);
    
    SF_A(i) = n * (W_mdd*l_rear_fuselage + (W_maa-W_mdd)*l_rear_fuselage/2 + W_tail_fin + W_engine);
    BM_A(i) = n*cos(deg2rad(result_alpha(i))-0.0159)*(l_rear_fuselage^2/2*W_mdd + l_rear_fuselage^2/6*(W_maa-W_mdd) + l_tail_fin*W_tail_fin * L_engine*W_engine); 
    % 0.0159 is the angle of incidence
    
    SF_C(i) = n*(W_mdd*(l_rear_fuselage-x_cc) + (W_mcc-W_mdd)*(l_rear_fuselage-x_cc)/2  + W_tail_fin + W_engine);
    BM_C(i) = n*cos(deg2rad(result_alpha(i))-0.0159)*((l_rear_fuselage-x_cc)^2/2*W_mdd + (l_rear_fuselage-x_cc)^2/6*(W_mcc-W_mdd) + (l_tail_fin-x_cc)*W_tail_fin * (l_engine-x_cc)*W_engine); 
    
    SF_B(i) = n*(W_mdd*(l_rear_fuselage-x_bb) + (W_mcc-W_mdd)*(l_rear_fuselage-x_cc)/2 + W_tail_fin + W_engine);
    BM_B(i) = n*cos(deg2rad(result_alpha(i))-0.0159)*((l_rear_fuselage-x_bb)^2/2*W_mdd + (l_rear_fuselage-x_bb)^2/6*(W_mbb-W_mdd)+ (l_tail_fin-x_bb)*W_tail_fin * (l_engine-x_bb)*W_engine);
    
    
    Ty_aa(i) = - F_fin;
    Tz_aa(i) = -(SF_A(i) - result_L_t(i))*cos(deg2rad(result_alpha(i))-0.0159);
    My_aa(i) =   BM_A(i) - result_L_t(i)*l_tail_fin*cos(deg2rad(result_alpha(i))-0.0159);
    Mz_aa(i) =   F_fin*l_tail_fin;
    Mx_aa(i) = - F_fin*abs(D.VT.y);
    
    
    Ty_cc(i) = -F_fin;
    Tz_cc(i) = -(SF_C(i) - result_L_t(i))*cos(deg2rad(result_alpha(i))-0.0159);
    My_cc(i) = BM_C(i) - result_L_t(i)*(l_tail_fin-x_cc)*cos(deg2rad(result_alpha(i))-0.0159);
    Mz_cc(i) = F_fin*(l_tail_fin-x_cc);
    Mx_cc(i) = - F_fin*abs(D.VT.y);
    
    
    Ty_bb(i) = -F_fin;
    Tz_bb(i) = -(SF_B(i) - result_L_t(i))*cos(deg2rad(result_alpha(i))-0.0159);
    My_bb(i) = BM_C(i) - result_L_t(i)*(l_tail_fin-x_bb)*cos(deg2rad(result_alpha(i))-0.0159);
    Mz_bb(i) = F_fin*(l_tail_fin-x_bb);
    Mx_bb(i) = - F_fin*abs(D.VT.y);

 end