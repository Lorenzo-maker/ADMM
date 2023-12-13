function data = fsae_multibody_vehicle_data_template(track, alfarange)

% Define the parameters of a multibody vehicle model of the FSAE type.

% Vehicle type of this dataset
data.model_type = 'FSAEMultibodyVehicle';

%% Suspension geometry

% The default numeric data below are read from a MSC Adams file describing a FSAE car. Numbers with 4 decimal digits are read from Adams file.
% These data are relative to a reference frame {D} attached to the differential, located at the midpoint of the rear axle.
% Numeric data are transformed to auxiliary frames; hub frames {H} are attached to the knuckles, and convenience frames {C} are attached to the chassis facing the wheels.
%
% The current suspension geometry uses a push-rod and rocker configuration. Other configurations are possible as well.
% The reference frame {R} is attached to the rocker.
%
% The chassis frame {B} is located at the center or mass of the chassis.
%
% All frames (except {R}) have the same orientation.
%
% All data are expressed in S.I. units

% Origin of chassis frame {B} w.r.t. differential {D} (fixed, but unused here)
% d_db = [ 0.815; 0.0000; 0.0100];

% Origin of interface frames {C} w.r.t. differential {D} (fixed)
d_dc_FL = [ 1.5800; 0.2000; 0.0000];
d_dc_RL = [ 0.0000; 0.2000; 0.0000];

% Origin of knuckle frames {H} w.r.t. differential {D} (at reference config in Adams file)
d_dh_FL = [ 1.5800; 0.6050; 0.0000];
d_dh_RL = [ 0.0000; 0.5550; 0.0000];

% Origin of rocker frames {R} w.r.t. differential {D} (at reference config in Adams file)
d_dr_FL = [ 1.5223; 0.2856; 0.0982];
d_dr_RL = [ 0.1050; 0.3019; 0.0433];

% Euler ZYX angles of rocker frames {R} w.r.t. differential {D} (at reference config in Adams file)
E_dr_FL = [ 1.1964;-0.2751; 1.7102];
E_dr_RL = [ 1.3882;-0.4945; 2.0827];

% Susp mounting points on chassis w.r.t. interface frames {C}
data.P_c_susp_FL = [
    [ 1.4195; 0.3016; 0.0289], ... % 1
    [ 1.7006; 0.2894; 0.0231], ... % 2
    [ 1.4186; 0.2522;-0.1293], ... % 3
    [ 1.6976; 0.2400;-0.1300], ... % 4
    [ 1.7400; 0.2450;-0.1103], ... % 5 (steering-rod)
    ]-d_dc_FL;
data.P_c_susp_RL = [
    [ 0.0388; 0.3085; 0.0273], ... % 1
    [ 0.1982; 0.3370; 0.0401], ... % 2
    [-0.0229; 0.2407;-0.1485], ... % 3
    [ 0.1929; 0.3232;-0.1378], ... % 4
    [-0.0553; 0.2214;-0.1508], ... % 5 (steering-rod)
    ]-d_dc_RL;

% Susp mounting points on knuckle w.r.t. knuckle frames {H}
data.P_h_susp_FL = [
    [ 1.5800; 0.5620; 0.1270], ... % 1
    [ 1.5800; 0.5620; 0.1270], ... % 2
    [ 1.5920; 0.5710;-0.1130], ... % 3
    [ 1.5920; 0.5710;-0.1130], ... % 4
    [ 1.6400; 0.5710;-0.0760], ... % 5 (steering-rod)
    ]-d_dh_FL;
data.P_h_susp_RL = [[ 0.0000; 0.5120; 0.1250], ... % 1
    [ 0.0000; 0.5120; 0.1250], ... % 2
    [ 0.0120; 0.5210;-0.1150], ... % 3
    [ 0.0120; 0.5210;-0.1150], ... % 4
    [-0.0230; 0.5340;-0.1135]] ... % 5 (steering-rod)
    -d_dh_RL;

% Spring mounting points on chassis w.r.t. interface frame {C}
data.P_c_spring_FL = [ 1.5107; 0.1748; 0.2879]-d_dc_FL;
data.P_c_spring_RL = [ 0.0361; 0.2184;-0.1141]-d_dc_RL;

% Spring mounting points on rocker w.r.t. rocker frame {R}
data.P_r_spring_FL = [ 0.0818; 0.1010; 0.0000];
data.P_r_spring_RL = [-0.0653; 0.0408; 0.0000];

% Push-rod mounting points on rocker w.r.t. rocker frame {R}
data.P_r_push_FL = [ 0.1540; 0.0000; 0.0000];
data.P_r_push_RL = [ 0.0680; 0.0000; 0.0000];

% Push-rod mounting points on knuckle w.r.t. knuckle frame {H}
data.P_h_push_FL = [ 1.5886; 0.5395;-0.0595]-d_dh_FL;
data.P_h_push_RL = [ 0.0320; 0.4901;-0.0783]-d_dh_RL;

% This defines the geometry of the suspensions.
% To modify the wheelbase/trackwidth/height of the car (while preserving the susp geometry), use the variables below.
%
% TODO We can create a GUI to enter these values.

% Adjust the distance of axles from CG
data.a1 = 0.765;
data.a2 = 0.815;

% Adjust the track widths (distance between the wheel centers)
data.t1 = 1.21;
data.t2 = 1.11;

% Adjust the height of CG above the axles (positive values shift the axles downward)
data.h1 = 0.01;
data.h2 = 0.01;

% Compute the axle widths (distance between interface frames)
data.d1 = data.t1-2*(d_dh_FL(2)-d_dc_FL(2));
data.d2 = data.t2-2*(d_dh_RL(2)-d_dc_RL(2));

% Compute the pose of interface frames {C} w.r.t. chassis frame {B}
data.q_bc_FL = [ data.a1;data.d1/2;-data.h1;0;0;0];
data.q_bc_RL = [-data.a2;data.d2/2;-data.h2;0;0;0];

% Compute the pose of knuckle frames {H} w.r.t. interface frames {C}
data.q_ch_FL = [d_dh_FL-d_dc_FL;0;0;0];
data.q_ch_RL = [d_dh_RL-d_dc_RL;0;0;0];

% Compute the pose of rocker frames {R} w.r.t. interface frames {C}
data.q_cr_FL = [d_dr_FL-d_dc_FL;E_dr_FL];
data.q_cr_RL = [d_dr_RL-d_dc_RL;E_dr_RL];

%% Mechanical properties of the susp

% Cell indices {i:1|2,j:1|2} refer to the {axle:F|R,side:L|R} wheel.
% For properties shared between wheels on the same axis/side, only the i/j index is used.

% Spring stiffness
data.k_s{1} = 35980;
data.k_s{2} = 24080;

% Damper damping
data.c_d{1} = 3280;
data.c_d{2} = 2195;

% Anti-roll bar stiffness
data.k_a{1} = 7750;
data.k_a{2} = 5000;

%% Steering geometry

% Direction of rack translation (when steering right)
data.steer_dir{1} = [0;1;0];
data.steer_dir{2} = [0;1;0];

% Steering ratio (rack translation over wheel rotation)
data.steer_ratio{1} = 0.0143;
data.steer_ratio{2} = 0; % NOTE Rear wheels do not steer.

% Maximum steering wheel rotation (half range)
data.steer_max = 2/3*pi;

%% Other properties of the susp

% Length of the springs at rest
data.coil_length{1} = 0.1952;
data.coil_length{2} = 0.2058;

% Maximum travel of the suspension trim
data.trim_lim{1} = [-0.05,0.05];
data.trim_lim{2} = [-0.05,0.05];

% Stiffness of the hard-stops
data.bump_order{1} = 11;
data.bump_order{2} = 11;

%% Inertial properties

% Mass of the wheels
data.m_w{1} = 10;
data.m_w{2} = 10;

% Inertia of the wheels
data.i_xx_w{1} = 0.4;
data.i_xx_w{2} = 0.4;
data.i_yy_w{1} = 0.7;
data.i_yy_w{2} = 0.7*1.5;

% Mass of the chassis
data.m_b = 241.38;

% Inertia of the chassis
data.i_xx_b = 41;
data.i_yy_b = 100;
data.i_zz_b = 110;
data.i_xz_b = 0;

%% Aerodynamics

% Vehicle frontal area
data.aero_S = 1.4;

% Air density
data.aero_rho = 1.2073;

% Shape coefficients
data.aero_Cx = 0.84;
data.aero_Cz = 1.34;

% Aero balance (front over total)
data.aero_balance = 0.4;

% Nominal height from ground (needed for moment of drag force according to Guiggiani)
data.h0 = 0.25;

%% Tire properties

% Load TIR data with tire parameters
load('fsae_2022.mat','tir')
data.tir{1} = tir;
data.tire_radius{1} = tir.UNLOADED_RADIUS;
data.tire_k{1} = tir.VERTICAL_STIFFNESS;
data.tire_c{1} = tir.VERTICAL_DAMPING;
load('fsae_2022.mat','tir')
data.tir{2} = tir;
data.tire_radius{2} = tir.UNLOADED_RADIUS;
data.tire_k{2} = tir.VERTICAL_STIFFNESS;
data.tire_c{2} = tir.VERTICAL_DAMPING;

%% Other properties of the vehicle

% Brake balance (front over total)
data.brake_balance = 0.6;

% Max driving power
data.Pa_max = 47000;

% Max braking power
data.Pb_max = 1800000;

% Max driving torque
data.Ta_max = 650;

% Max braking torque
data.Tb_max = 2500;

% Gravity
data.g = 9.81;


%% Compute track quantities (necessary to define scaling factors)

import casadi.*

grid = alfarange;
N = length(grid)-1;

% Extract data
tw_l = zeros(1,N+1);
tw_r = zeros(1,N+1);
d_gs = zeros(3,N+1);
R_gs = zeros(9,N+1);
pos_map = track.fun_pos.map(N+1, 'thread', 4);
pos = pos_map(grid);
d_gs = full(pos(1:3,:));
tw_l = pos(4,:);
tw_r = pos(5,:);
alpha_sym = SX.sym('alpha_sym',1);
Rgs_vec = Function('Rgs_vec', {alpha_sym}, {reshape(full(track.fun_Rgs(alpha_sym)),9,1)});
Rgs_vec_map = Rgs_vec.map(N+1, 'thread', 4);
R_gs = full(Rgs_vec_map(grid));
% for k = 1:N+1
%     tmp = full(track.fun_pos(grid(k)));
%     d_gs(:,k) = tmp(1:3);
%     tw_l(k) = tmp(4);
%     tw_r(k) = tmp(5);
%     tmp = full(track.fun_Rgs(grid(k)));
%     R_gs(:,k) = reshape(tmp,9,1);
% end


% Compute length of the specified track segment
len = sum(sqrt(sum((d_gs(:,2:N+1)-d_gs(:,1:N)).^2)));

% Save data to output structure
data.track = struct;
data.track.grid = grid;
data.track.len = len;
data.track.d_gsR_gs = [d_gs;R_gs];
data.track.width = [tw_l;tw_r];
data.track.psi_grid = atan_track_diff(R_gs(1:3,:));

lap = all((diff(alfarange)<0)) + 1;



%% Scale factors

% Scale factors for states
data.X_scale = [
    50;
    10;
    5;
    1;
    1;
    2;
    1.5;
    1.5;
    1.5;
    1.5;
    200;
    200;
    200;
    200;
    max(abs(d_gs(1,:)))+100;%5000; % NOTE It can be refined after track is loaded.
    max(abs(d_gs(2,:)))+100; % NOTE It can be refined after track is loaded.
    max(abs(d_gs(3,:)))+100; % NOTE It can be refined after track is loaded.
    (lap+1)*2*pi; % NOTE It can be refined after track is loaded.
    pi/4; % NOTE It can be refined after track is loaded.
    pi/4; % NOTE It can be refined after track is loaded.
    0.1;
    0.1;
    0.1;
    0.1;
    ];

% Scale factors for inputs
data.U_scale = [
    data.Ta_max;
    data.Tb_max;
    data.steer_max;
    ];

% Scale factors for algebraic variables
F_tire_scale = [4000,4000,4000];
data.Z_scale = [
    10; % NOTE It can be refined after track is loaded.
    1;
    1;
    reshape(repmat(F_tire_scale,4,1),[],1);
    ];

end
