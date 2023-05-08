 clear all; close all; clc;
addpath(genpath('..\..\Casadi'))
addpath('C:\Users\blore\OneDrive - University of Pisa\Unipi\PhD\NURBS\Track')
addpath('C:\Users\blore\OneDrive - University of Pisa\Unipi\PhD\Tyre\MF_Tyre')
addpath('..\..\Functions_utils')
addpath('..\')
addpath('..\..\')
import casadi.*

%% Track & Grid
warning off
load('track_stored.mat')
warning on
track = track_stored{2};
% track = track_fun(1, '2D', 1, 'end', 200); 
track.residual
track.tot_error

% Grid parameters
N = 600;
d = 2;

% Computed grid
%alfa_grid = linspace(0.55,0.85,N+1);
alfa_grid = linspace(0.,0.1,N+1);

dalfa = alfa_grid(2) - alfa_grid(1);
% alfa in collocation points
import casadi.*
tau_root = [0 collocation_points(d, 'legendre')];
alfa_grid_colloc = zeros(d+1, N+1);
alfa_grid_colloc(1,:) = alfa_grid';
for i = 2:length(tau_root)
    alfa_grid_colloc(i,:) = alfa_grid + dalfa*tau_root(i);
end
pos_grid = full(track.fun_pos(alfa_grid));
ns_grid  = full(track.fun_vh(alfa_grid));
kp_grid  = full(track.fun_kp(alfa_grid));
ts_grid = full(track.fun_ts(alfa_grid));
ts_grid(3,:) = 0;
ts_grid = ts_grid./vecnorm(ts_grid);
psi_grid = atan_track(ts_grid, 'clockwise');
%
psi_grid_colloc = zeros(d, N+1);
for i = 1:d
    ts_grid_colloc = full(track.fun_ts(alfa_grid_colloc(i+1,:)));
    ts_grid_colloc(3,:) = 0;
    ts_grid_colloc = ts_grid_colloc./vecnorm(ts_grid_colloc);
    psi_grid_colloc(i,:) = atan_track(ts_grid_colloc, 'clockwise');
end
%
RVSG = cell(1,N+1);
RVSG_colloc = cell(d,N+1);
for i = 1:N+1
    RVSG{i} = [cos(-psi_grid(i)), sin(-psi_grid(i)), 0; -sin(-psi_grid(i)), cos(-psi_grid(i)), 0; 0, 0, 1]*(full(track.fun_Rgs(alfa_grid(i))))';
    for j = 1:d
        RVSG_colloc{j, i} = [cos(-psi_grid_colloc(j,i)), sin(-psi_grid_colloc(j,i)), 0; -sin(-psi_grid_colloc(j,i)), cos(-psi_grid_colloc(j,i)), 0; 0, 0, 1]*(full(track.fun_Rgs(alfa_grid_colloc(j+1,i))))';
    end
end
RVSG_colloc = RVSG_colloc(:, 1:end-1); % N dimension

%% Tyre
load('Fsae_2022.mat')
[Fx0, Fy0, data.tyre.Fx, data.tyre.Fy, data.tyre.Fxlim, data.tyre.Fylim, data.tyre.Gxa, data.tyre.Gyk] = MF_Tyre(tir, 'shift', 0, 'curvature', 0);

% vehicle data
data.rho         = 1.2073;                    % air density
data.S           = 1.4;                       % front area
data.cx          = 0.84;                      % drag coefficient
data.cz          = 1.34;                      % lift coefficient
data.ab          = 0.4;
data.a1          = 0.765;
data.a2          = 0.815;
data.l           = data.a1 + data.a2;
data.hg          = 0.3;%0.28;
data.t1          = 1.2;
data.t2          = 1.2;
data.rw          = tir.UNLOADED_RADIUS;
data.m           = 240;                       % mass of car in kg
data.Jw          = (1+3)*0.5*tir.MASS*tir.UNLOADED_RADIUS^2;
data.Jz          = data.m*data.a1*data.a2*0.92;
data.G           = 9.81;                      % gravity
data.Fz_start    = data.m*data.G;             % weigth force (Newton)
data.mu_x        = 3.0 ;                      % dry road example
data.mu_y        = 3.0;                       % dry road example
data.tol         = 10^-3;                     % tolerance for avoiding NaN
data.Vm          = 22.23;                     % mean total speed (m/s) ~ 80 km/h
data.hm          = 0.3;                       % mean discretization 
data.Vi          = 30;                        % initial speed
data.Vmax        = 150/3.6;                   % max speed
data.Pmin        = -1800e3;                   % max braking power
data.Pmax        = 47e3;                      % max power
data.kb          = 0.6;
data.q1          = (2/5)*tir.UNLOADED_RADIUS;%0.06;
data.q2          = (1/5)*tir.UNLOADED_RADIUS;%0.065;
data.qm          = (data.a2*data.q1 + data.a1*data.q2)/data.l;
data.k_phi1      = 0.5*36000*data.t1^2;
data.k_phi2      = 0.5*24000*data.t2^2;
data.k_phi       = (data.k_phi1 + data.k_phi2 + 1e-5);
data.eta         = 1;
data.tau = 1/3;
data.d0 = 0;
data.ep1 = 0;

% Scaling Factors
data.u_scale  = data.Vmax;
data.v_scale  = 10;
data.r_scale  = 8;
data.x_max = max(pos_grid(1,:))+100;
data.x_min = min(pos_grid(1,:))-100;
data.x_scale = max(abs(pos_grid(1,:)));
data.y_max = max(pos_grid(2,:))+100;
data.y_min = min(pos_grid(2,:))-100;
data.y_scale = max(abs(pos_grid(2,:)));
data.psi_scale = 4*pi;
data.w_scale = data.Vmax/data.rw;
data.Ta_scale = 2000;%data.Pmax/(data.w_scale);
data.Tb_scale = 10000;%data.Pmin/(data.w_scale);
data.delta_scale = pi;
data.Fz_scale = data.m*data.G + 0.5*data.rho*data.S*data.cz*data.Vmax^2;
data.Fx_scale = data.Fz_scale*data.mu_x;
data.Fy_scale = data.Fz_scale*data.mu_y;
data.P_scale  = max(abs(data.Pmin), data.Pmax);
data.ep_scale = 10;
data.h_scale = 1;
data.X_scale = [data.u_scale; data.v_scale; data.r_scale; data.x_scale; data.y_scale; data.psi_scale; repmat(data.w_scale,4,1)];
data.U_scale = [data.Ta_scale; data.Tb_scale; data.delta_scale];
data.Z_scale = [repmat(data.Fz_scale,4,1);repmat(data.Fx_scale,4,1);repmat(data.Fy_scale,4,1); data.ep_scale; data.h_scale];

% Car Model
car = vehicle_casadi('double-track-full', data);
% Intial Condition
if pos_grid(1,1) == 0
    psi_0 = -pi;
    y0_lb = -pos_grid(4,1);
    y0_ub =  pos_grid(4,1);
else
    psi_0 = psi_grid(1);
    y0_lb = pos_grid(2,1);
    y0_ub = pos_grid(2,1);
end
car.Xb.x0 = [car.data.Vi;0;0;pos_grid(1,1);pos_grid(2,1);psi_grid(1);repmat(car.data.Vi/car.data.rw,2,1); repmat(car.data.Vi*(1 - 0)/car.data.rw,2,1)]./car.data.X_scale;
car.Xb.x0_lb = [car.data.Vi;0;0;pos_grid(1,1);y0_lb;psi_0;repmat(car.data.Vi/car.data.rw,2,1); repmat(car.data.Vi*(1 - 0.05)/car.data.rw,2,1)]./car.data.X_scale;
car.Xb.x0_ub = [car.data.Vi;0;0;pos_grid(1,1);y0_ub;psi_0;repmat(car.data.Vi/car.data.rw,2,1); repmat(car.data.Vi*(1 + 0.05)/car.data.rw,2,1)]./car.data.X_scale;
% Numerical data for variable constraints
Pdata.Pgb = [-pos_grid(4,2:end); pos_grid(4,2:end)];
Rvsg = RVSG;
Rvsg = horzcat(Rvsg{:});
Rvsg = Rvsg(:);
Rvsg = reshape(Rvsg, 9, []);
Rvsg_colloc = horzcat(RVSG_colloc{:});
Rvsg_colloc = Rvsg_colloc(:);
Rvsg_colloc = reshape(Rvsg_colloc, 9*d, []);
Pdata.Pg = [pos_grid(1:2,2:end); ns_grid(1:2, 2:end); Rvsg(:,2:end); Rvsg_colloc];

% NLP
pb = nlp_casadi_2(car.nx, car.nu, car.nz, car.nuc, car.nzc, car.nwg, N, d, car.Xb, car.Ub, car.Zb, [], [], [], 4 + 9 + 9*d, 2, Pdata, 1, 0, 'SXMX');
pb.w0 = [];
for k = 1:N
    initial_guess_k;    
    pb.w0 = [pb.w0;w0_k];
end
pb.w0 = [car.Xb.x0; pb.w0];

% NLP
pb.append_g(pb.x_1(9) - pb.x_1(10), 0, 0, 'start');
RvsgKj = reshape(pb.pg(5:end),9, d+1);
RvsgK = reshape(RvsgKj(:, 1),3,3);
for i = 1:d   
   Rvsgj = reshape(RvsgKj(:, i+1),3,3);
   g3DKj = car.fun_g3D(pb.xc(6,i), Rvsgj); 
   pb.append_g(car.F(pb.xc(:,i), pb.u, pb.z, g3DKj)*pb.z(end)*car.data.Z_scale(end) - pb.xp(:,i).*car.data.X_scale, zeros(car.nx,1), zeros(car.nx,1)) ;
end
g3DK = car.fun_g3D(pb.x(6), RvsgK); 
xdot = car.F(pb.x, pb.u, pb.z, g3DK);
pb.Jj = (pb.z(end)*car.data.Z_scale(end))^2 + 10*((pb.u(3) - pb.u_1(3))*car.data.tau)^2 + 0.0001*xdot(2)^2 + 0.01*xdot(3)^2 + 0*0.0000001*(xdot(end-4)^2 + xdot(end-3)^2 + xdot(end-2)^2 + xdot(end-1)^2);
% Continuity equation for generic interval
pb.append_g(pb.xc_end - pb.x , zeros(car.nx,1), zeros(car.nx,1));
% Stay on track constraint
pb.append_g(pb.x(4:5) - (pb.pg(1:2) + pb.pg(3:4)*pb.z(end-1)*car.data.Z_scale(end-1))./car.data.X_scale(4:5), zeros(2,1), zeros(2,1));
pb.append_g(pb.z(end-1), pb.pgb(1)/car.data.Z_scale(end-1), pb.pgb(2)/car.data.Z_scale(end-1));    
% Algebraic equations constraints
pb.append_g(car.Eq(pb.x, pb.u, pb.z, g3DK), zeros(car.neq,1), zeros(car.neq,1));
% Steering constraint
% pb.append_g(pb.u(3) - pb.u_1(3), -pi/data.U_scale(3), pi/data.U_scale(3));
% Complementary constraint
pb.append_g(pb.u(1)*pb.u(2), -1e-6, 1e-6);
% PowerTrain
pb.append_g(((pb.u(1)*car.data.U_scale(1)*pb.x(1)*car.data.X_scale(1)/tir.UNLOADED_RADIUS) + (pb.u(2)*car.data.U_scale(2)*pb.x(1)*car.data.X_scale(1)/tir.UNLOADED_RADIUS))/car.data.P_scale, car.data.Pmin/car.data.P_scale, car.data.Pmax/car.data.P_scale);
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


pb.build_solver(opts, 'NLP');
pb.solution

% Extract solution
U_opt = pb.sol_num(1:pb.nu, :).*car.data.U_scale;
Z_opt = pb.sol_num(pb.nu + 1: pb.nu + pb.nz, :).*car.data.Z_scale;
X_opt = [pb.X0_num, pb.sol_num(pb.nu + pb.nz + pb.nx*d + 1:end, :)].*car.data.X_scale;
Qn = full(car.Q(X_opt(:, 2:end)./data.X_scale, U_opt./data.U_scale, Z_opt./data.Z_scale));
Fx11_opt = full(data.tyre.Fx(Qn(1,:), Qn(5,:), Z_opt(1,:), 0));
Fx12_opt = full(data.tyre.Fx(Qn(2,:), Qn(6,:), Z_opt(2,:), 0));
Fx21_opt = full(data.tyre.Fx(Qn(3,:), Qn(7,:), Z_opt(3,:), 0));
Fx22_opt = full(data.tyre.Fx(Qn(4,:), Qn(8,:), Z_opt(4,:), 0));

Fy11_opt = full(data.tyre.Fy(Qn(5,:), Qn(1,:), Z_opt(1,:), 0));
Fy12_opt = full(data.tyre.Fy(Qn(6,:), Qn(2,:), Z_opt(2,:), 0));
Fy21_opt = full(data.tyre.Fy(Qn(7,:), Qn(3,:), Z_opt(3,:), 0));
Fy22_opt = full(data.tyre.Fy(Qn(8,:), Qn(4,:), Z_opt(4,:), 0));

d11 = -data.d0 + data.tau*U_opt(3,:) + data.ep1*data.t1*(data.tau*U_opt(3,:)).^2/(2*data.l);
d12 =  data.d0 + data.tau*U_opt(3,:) - data.ep1*data.t1*(data.tau*U_opt(3,:)).^2/(2*data.l);

sol.U = U_opt;
sol.X = X_opt;
sol.Z = Z_opt;
sol.Fx = [Fx11_opt; Fx12_opt; Fx21_opt; Fx22_opt];
sol.Fy = [Fy11_opt; Fy12_opt; Fy21_opt; Fy22_opt];

g3D_opt = zeros(3,N+1);
for i = 1:length(X_opt(6,:))
    g3D_opt(:, i) = full(car.fun_g3D(full(X_opt(6,i)),RVSG{i})); 
end


% Optimal time
time_opt = [0, cumsum(Z_opt(end,:))];
disp(['Optimal time is: ', num2str(time_opt(end))])

close all
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
% Trajectory plot
track.full_plot(alfa_grid, 'patch', false, 'dotted', false);
hold on
grid on
plot(X_opt(4,:), X_opt(5,:), 'Linewidth', 2, 'color', [0, .8, .1]);


% Speed
figure(2)
plot(alfa_grid, X_opt(1,:))
title('u')

% wheel
figure(3)
plot(alfa_grid, X_opt(7,:), alfa_grid, X_opt(8,:), 'o', alfa_grid, X_opt(9,:), alfa_grid, X_opt(10,:),'*')
title('wij')
legend('$w_11$','$w_12$','$w_21$','$w_22$','Orientation','horizontal')

% dt
figure(4)
plot(alfa_grid(1:end-1), Z_opt(end,:))
title('h')

% Ta-Tb
figure(5)
plot(alfa_grid(1:end-1), U_opt(1,:), alfa_grid(1:end-1), U_opt(2,:))
title('T')

figure(7)
plot(alfa_grid(1:end-1), Qn(1,:))
hold on
plot(alfa_grid(1:end-1), Qn(2,:), 'o')
plot(alfa_grid(1:end-1), Qn(3,:))
plot(alfa_grid(1:end-1), Qn(4,:), '*')
title('k-slip')
legend

figure(8)
plot(X_opt(2,:))
hold on
plot(X_opt(3,:), 'o')
title('v-r')
legend

figure(9)
plot(alfa_grid(1:end-1), U_opt(3,:))
title('delta')

figure(10)
plot(alfa_grid(1:end-1), Z_opt(1,:), alfa_grid(1:end-1), Z_opt(2,:), 'o', alfa_grid(1:end-1), Z_opt(3,:), alfa_grid(1:end-1), Z_opt(4,:), '*')
title('Fzij')

figure(11)
plot(alfa_grid(1:end-1), Fx11_opt, alfa_grid(1:end-1), Fx12_opt, 'o', alfa_grid(1:end-1), Fx21_opt, alfa_grid(1:end-1), Fx22_opt, '*')
title('Fxij')
legend('$Fx_11$','$Fx_12$','$Fx_21$','$Fx_22$','Orientation','horizontal')
% 
% figure(12)
% plot(Fx11_opt, Fy11_opt, '*', Fx12_opt, Fy12_opt, 'o')


% [Fxmin21, Fxmax21] = (data.tyre.Fxlim(Z_opt(3,:),0));
% Fxmin21 = data.tyre.Gxa(Qn(3,:), Qn(7,:), Z_opt(3,:), 0).*Fxmin21;
% Fxmax21 = data.tyre.Gxa(Qn(3,:), Qn(7,:), Z_opt(3,:), 0).*Fxmax21;
% 
% [Fymin11, Fymax11] = (data.tyre.Fylim(Z_opt(1,:),0));
% Fymin11 = data.tyre.Gxa(Qn(5,:), Qn(1,:), Z_opt(1,:), 0).*Fymin11;
% Fymax11 = data.tyre.Gxa(Qn(5,:), Qn(1,:), Z_opt(1,:), 0).*Fymax11;
% Fxnorm11 = zeros(length(Fx11_opt),1);
% Fynorm11 = zeros(length(Fx11_opt),1);
% for i = 1:length(Fx11_opt)
%     Fxnorm11(i) = full(if_else_smooth(Fx11_opt(i), 0, Fx11_opt(i)./Fxmax11(i),  -Fx11_opt(i)./Fxmin11(i)));
%     Fynorm11(i) = full(if_else_smooth(Fy11_opt(i), 0, Fy11_opt(i)./Fymax11(i),  -Fy11_opt(i)./Fymin11(i)));
% end

%%
kmin = tir.KPUMIN;
kmax = tir.KPUMAX;
al_min = tir.ALPMIN;
al_max =  tir.ALPMAX;

figure(13);
ax1 = subplot(4,4,3);
axis(ax1, 'equal')
ax2 = subplot(4,4,4);
ax3 = subplot(4,4,7);
ax4 = subplot(4,4,8);
axS = subplot(4,4,[1,2,5,6]);
title(axS, 'Speed')
hold(axS, 'on')   
grid(axS,'on')
plot(alfa_grid, X_opt(1,:), 'parent', axS)  
point_speed = plot(axS, nan, nan, 'marker', 'o', 'color', 'b','HandleVisibility','off');
axT = subplot(4,4,9:16);
title(axT, 'Trajectory')
hold(axT, 'on')   
grid(axT,'on')
track.full_plot(alfa_grid, 'patch', false, 'dotted', false, 'Parent', axT)
plot(X_opt(4,:), X_opt(5,:), 'color', [0 0.5 0.5], 'Linewidth', 2);
tr = hgtransform('parent', axT);
pl = plot(0,0, 'parent', tr, 'marker', 'o', 'color', 'b');

al_grid1 = nonlinspace(0, al_max, 21, 0.6);
al_grid1 = [fliplr(-al_grid1), al_grid1];
k_grid1 = linspace(kmin,kmax,11);
al_grid2 = linspace(al_min, al_max, 11);
k_grid2 = nonlinspace(0, kmax, 21, 0.6);
k_grid2 = [fliplr(-k_grid2), k_grid2];

% n = 1;
% for j = 250:2:length(Fx11_opt)
%         [Fxmin, Fxmax] = data.tyre.Fxlim(Z_opt(1,j),0);
%         FxmaxM(n) = Fxmax;
%         FxmaxM = [FxmaxM; Fxmax];
%         FymaxM(n) = Fxmax;
%         [Fymin, Fymax] = data.tyre.Fylim(Z_opt(1,j),0);
%         plot(full(data.tyre.Fx(k_grid1(i), al_grid1, Z_opt(1,j), 0)/Fxmax),  full(data.tyre.Fy(al_grid1, k_grid1(i), Z_opt(1,j), 0)/Fymax), 'Parent', ax1)                

for j = 1:5:length(Fx11_opt)
    %%%% Wheel 11
    hold(ax1, 'on')   
    grid(ax1,'on')
    title(ax1, '$w_{11}$')
    [Fxmin, Fxmax] = data.tyre.Fxlim(Z_opt(1,j),0);
    [Fymin, Fymax] = data.tyre.Fylim(Z_opt(1,j),0);    
    for i = 1:length(k_grid1)
        plot(full(data.tyre.Fx(k_grid1(i), al_grid1, Z_opt(1,j), 0)/Fxmax),  full(data.tyre.Fy(al_grid1, k_grid1(i), Z_opt(1,j), 0)/Fymax), 'Parent', ax1)                
    end
    for i = 1:length(al_grid2)
        plot(full(data.tyre.Fx(k_grid2, al_grid2(i), Z_opt(1,j), 0)/Fxmax),  full(data.tyre.Fy(al_grid2(i), k_grid2, Z_opt(1,j), 0)/Fymax), 'Parent', ax1)
    end
    if j > 10
        [~,Fxmax10] = data.tyre.Fxlim(Z_opt(1,j-10: j),0);
        [~,Fymax10] = data.tyre.Fylim(Z_opt(1,j-10: j),0);       
        plot(Fx11_opt(j-10:j)./full(Fxmax10), Fy11_opt(j-10:j)./full(Fymax10), 'g-o', 'Parent', ax1)
    else
        plot(Fx11_opt(j)/full(Fxmax), Fy11_opt(j)/full(Fymax), 'go', 'Parent', ax1)
    end
    %%%%% Wheel 12
    hold(ax2, 'on')
    grid(ax2, 'on')
    title(ax2, '$w_{12}$')
    [Fxmin, Fxmax] = data.tyre.Fxlim(Z_opt(2,j),0);
    [Fymin, Fymax] = data.tyre.Fylim(Z_opt(2,j),0);    
    for i = 1:length(k_grid1)
        plot(full(data.tyre.Fx(k_grid1(i), al_grid1, Z_opt(2,j), 0)/Fxmax),  full(data.tyre.Fy(al_grid1, k_grid1(i), Z_opt(2,j), 0)/Fymax), 'Parent', ax2)                
    end
    for i = 1:length(al_grid2)
        plot(full(data.tyre.Fx(k_grid2, al_grid2(i), Z_opt(2,j), 0)/Fxmax),  full(data.tyre.Fy(al_grid2(i), k_grid2, Z_opt(2,j), 0)/Fymax), 'Parent', ax2)
    end
    if j > 10
        [~,Fxmax10] = data.tyre.Fxlim(Z_opt(2,j-10: j),0);
        [~,Fymax10] = data.tyre.Fylim(Z_opt(2,j-10: j),0);       
        plot(Fx12_opt(j-10:j)./full(Fxmax10), Fy12_opt(j-10:j)./full(Fymax10), 'g-o', 'Parent', ax2)
    else
        plot(Fx12_opt(j)/full(Fxmax), Fy12_opt(j)/full(Fymax), 'go', 'Parent', ax2)
    end
    %%%%% Wheel 21
    hold(ax3, 'on')   
    grid(ax3, 'on')
    title(ax3, '$w_{21}$')   
    [Fxmin, Fxmax] = data.tyre.Fxlim(Z_opt(3,j),0);
    [Fymin, Fymax] = data.tyre.Fylim(Z_opt(3,j),0);    
    for i = 1:length(k_grid1)
        plot(full(data.tyre.Fx(k_grid1(i), al_grid1, Z_opt(3,j), 0)/Fxmax),  full(data.tyre.Fy(al_grid1, k_grid1(i), Z_opt(3,j), 0)/Fymax), 'Parent', ax3)                
    end
    for i = 1:length(al_grid2)
        plot(full(data.tyre.Fx(k_grid2, al_grid2(i), Z_opt(3,j), 0)/Fxmax),  full(data.tyre.Fy(al_grid2(i), k_grid2, Z_opt(3,j), 0)/Fymax), 'Parent', ax3)
    end
    if j > 10
        [~,Fxmax10] = data.tyre.Fxlim(Z_opt(3,j-10: j),0);
        [~,Fymax10] = data.tyre.Fylim(Z_opt(3,j-10: j),0);       
        plot(Fx21_opt(j-10:j)./full(Fxmax10), Fy21_opt(j-10:j)./full(Fymax10), 'g-o', 'Parent', ax3)
    else
        plot(Fx21_opt(j)/full(Fxmax), Fy21_opt(j)/full(Fymax), 'go', 'Parent', ax3)
    end    
    %%%%% Wheel 22
    hold(ax4,'on')
    grid(ax4, 'on')
    title(ax4, '$w_{22}$')  
    [Fxmin, Fxmax] = data.tyre.Fxlim(Z_opt(4,j),0);
    [Fymin, Fymax] = data.tyre.Fylim(Z_opt(4,j),0);    
    for i = 1:length(k_grid1)
        plot(full(data.tyre.Fx(k_grid1(i), al_grid1, Z_opt(4,j), 0)/Fxmax),  full(data.tyre.Fy(al_grid1, k_grid1(i), Z_opt(4,j), 0)/Fymax), 'Parent', ax4)                
    end
    for i = 1:length(al_grid2)
        plot(full(data.tyre.Fx(k_grid2, al_grid2(i), Z_opt(4,j), 0)/Fxmax),  full(data.tyre.Fy(al_grid2(i), k_grid2, Z_opt(4,j), 0)/Fymax), 'Parent', ax4)
    end
    if j > 10
        [~,Fxmax10] = data.tyre.Fxlim(Z_opt(4,j-10: j),0);
        [~,Fymax10] = data.tyre.Fylim(Z_opt(4,j-10: j),0);       
        plot(Fx22_opt(j-10:j)./full(Fxmax10), Fy22_opt(j-10:j)./full(Fymax10), 'g-o', 'Parent', ax4)
    else
        plot(Fx22_opt(j)/full(Fxmax), Fy22_opt(j)/full(Fymax), 'go', 'Parent', ax4)
    end  
    %%%
    hold(axS, 'on')
    grid(axS, 'on')
    set(point_speed, 'XData', alfa_grid(j), 'YData', X_opt(1,j))
    %%%
    hold(axT, 'on')
    grid(axT, 'on')
    view(axT, 0, 90)    
    tr.Matrix = diag(ones(4,1));
    tr.Matrix(1:3,4) = [X_opt(4,j), X_opt(5,j), 0];
    axis('equal')
    drawnow
    cla(ax1, 'reset')
    cla(ax2, 'reset')
    cla(ax3, 'reset')
    cla(ax4, 'reset')
end
