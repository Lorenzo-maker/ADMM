clear all; close all; clc;
addpath(genpath('..\..\Casadi'))
addpath('C:\Users\blore\OneDrive - University of Pisa\Unipi\PhD\NURBS\Track')
addpath('..\..\Functions_utils')
addpath('..\')
addpath('..\..\')
import casadi.*

% Load track
warning off
load('track_stored.mat')
warning on
% track = track_stored{1};
track = track_fun(1, '2D', 1, 'end', 200); 
track.residual
track.tot_error

% Grid parameters
N = 500;
d = 2;

% Computed grid
alfa_grid = linspace(0,1,N+1);
pos_grid = full(track.fun_pos(alfa_grid));
ns_grid  = full(track.fun_vh(alfa_grid));
kp_grid  = full(track.fun_kp(alfa_grid));
ts_grid = full(track.fun_ts(alfa_grid));
ts_grid(3,:) = 0;
ts_grid = ts_grid./vecnorm(ts_grid);
psi_grid = atan_track(ts_grid, 'clockwise');

% vehicle data
data.rho         = 1.2073;                    % air density
data.S           = 1.4;                       % front area
data.cx          = 1.19;                      % drag coefficient
data.cz          = 3.63;                      % lift coefficient
data.G           = 9.81;                      % gravity
data.m           = 647.7;                     % mass of car in kg
data.Fz_start    = data.m*data.G;                       % weigth force (Newton)
data.mu_x        = 1.0 ;                      % dry road example
data.mu_y        = 1.0;                       % dry road example
data.tol         = 10^-3;                     % tolerance for avoiding NaN
data.Vm          = 22.23;                     % mean total speed (m/s) ~ 80 km/h
data.hm          = 0.25;                      % mean discretization 
data.Vi          = 10;                        % initial speed
data.Vmax        = 320/3.6;                   % max speed
data.Pmin        = -1800e3;                   % max braking power
data.Pmax        = 470e3;                     % max power

% Scaling Factors
data.u_scale  = data.Vmax;
data.v_scale  = 6;
data.r_scale  = 4;
data.x_max = max(pos_grid(1,:))+100;
data.x_min = min(pos_grid(1,:))-100;
data.x_scale = max(abs(pos_grid(1,:)));
data.y_max = max(pos_grid(2,:))+100;
data.y_min = min(pos_grid(2,:))-100;
data.y_scale = max(abs(pos_grid(2,:)));
data.psi_scale = 4*pi;
data.Fz_scale = data.m*data.G + 0.5*data.rho*data.S*data.cz*data.Vmax^2;
data.Fx_scale = data.Fz_scale*data.mu_x;
data.rp_scale = 10;
data.Fy_scale = data.Fz_scale*data.mu_y;
data.P_scale  = max(abs(data.Pmin), data.Pmax);
data.ep_scale = 10;
data.h_scale = 5;
data.X_scale = [data.u_scale; data.v_scale; data.r_scale; data.x_scale; data.y_scale; data.psi_scale];
data.U_scale = [data.Fx_scale; data.rp_scale];
data.Z_scale = [data.Fz_scale; data.Fy_scale; data.P_scale; data.ep_scale; data.h_scale];

% Car Model
car = vehicle_casadi('point-mass', data);
% Intial Condition
car.Xb.x0 = [car.data.Vi;0;0;pos_grid(1,1);pos_grid(2,1);psi_grid(1)]./car.data.X_scale;
car.Xb.x0_lb = [car.data.Vi;0;0;pos_grid(1,1);-pos_grid(4,1);psi_grid(1)]./car.data.X_scale;
car.Xb.x0_ub = [car.data.Vi;0;0;pos_grid(1,1);pos_grid(4,1);psi_grid(1)]./car.data.X_scale;
% Numerical data for variable constraints
Pdata.Pgb = [-pos_grid(4,2:end); pos_grid(4,2:end)];
Pdata.Pg = [pos_grid(1:2,2:end); ns_grid(1:2, 2:end)];

% NLP
pb = nlp_casadi_2(car.nx, car.nu, car.nz, car.nuc, car.nzc, car.nwg, N, d, car.Xb, car.Ub, car.Zb, [], [], [], 4, 2, Pdata, 1, 0, 'SX');
pb.w0 = [];
for k = 1:N
    initial_guess_k;    
    pb.w0 = [pb.w0;w0_k];
end
pb.w0 = [car.Xb.x0; pb.w0];

% NLP
for i = 1:d
   pb.append_g(car.F(pb.xc(:,i), pb.u, pb.z)*pb.z(5)*car.data.Z_scale(5) - pb.xp(:,i).*car.data.X_scale, zeros(car.nx,1), zeros(car.nx,1)) ;
end
xdot = car.F(pb.x, pb.u, pb.z);
pb.Jj = (pb.z(5)*car.data.Z_scale(5))^2 + 0.1*(pb.u(2)*car.data.U_scale(2))^2 + 0.01*xdot(2)^2 + 1*(pb.z_1(2)-pb.z(2))^2;
% Evaluate cost function for generic interval
pb.cost;
% Continuity equation for generic interval
pb.append_g(pb.xc_end - pb.x , zeros(car.nx,1), zeros(car.nx,1));
% pb.append_g((pb.x - pb.x_1).*data.X_scale - xdot.*pb.z(5)*data.Z_scale(5), zeros(car.nx,1), zeros(car.nx,1))
% Stay on track constraint
pb.append_g(pb.x(4:5) - (pb.pg(1:2) + pb.pg(3:4)*pb.z(4)*car.data.Z_scale(4))./car.data.X_scale(4:5), zeros(2,1), zeros(2,1));
pb.append_g(pb.z(4), pb.pgb(1)/car.data.Z_scale(4), pb.pgb(2)/car.data.Z_scale(4));    
% Adherence constraint
pb.append_g((pb.z(2)*car.data.Z_scale(2)/(car.data.mu_y*pb.z(1)*car.data.Z_scale(1)))^2 + (pb.u(1)*car.data.U_scale(1)/(car.data.mu_x*pb.z(1)*car.data.Z_scale(1)))^2,0,1); 
% Algebraic equations constraints
pb.append_g(car.Eq(pb.x, pb.u, pb.z), zeros(car.neq,1), zeros(car.neq,1));
% Build map function 
pb.build_map;
% Build g 
pb.build_g;
% Build J
pb.build_J;


opts = struct;
opts.ipopt.linear_solver = 'ma57';
opts.ipopt.mu_init = 1e-3;
% opts.ipopt.mu_strategy = 'adaptive';
opts.ipopt.max_iter = 5000;
opts.ipopt.print_level = 5;


pb.build_solver([], 'NLP');
pb.solution

% Extract solution
U_opt = pb.sol_num(1:pb.nu, :).*car.data.U_scale;
Z_opt = pb.sol_num(pb.nu + 1: pb.nu + pb.nz, :).*car.data.Z_scale;
X_opt = [pb.X0_num, pb.sol_num(pb.nu + pb.nz + pb.nx*d + 1:end, :)].*car.data.X_scale;

% Optimal time
time_opt = [0, cumsum(Z_opt(end,:))];
disp(['Optimal time is: ', num2str(time_opt(end))])

% Trajectory plot
track.full_plot(alfa_grid, 'patch', false, 'dotted', false)
hold on
grid on
plot(X_opt(4,:), X_opt(5,:), 'Linewidth', 2, 'color', [0, .8, .1])

% Speed
figure(2)
plot(alfa_grid, X_opt(1,:))
