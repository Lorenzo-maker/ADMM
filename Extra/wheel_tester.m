% script simulazione slittamento ruota
close all; clear all;clc;
addpath('Functions_utils')
% T = 3; % torque (for constant torque simulation)
r = 0.25; % wheel radius
m = 30;   % wheel mass

% Magic formula coefficients
B = 20;
C = 1.5;
D = 2;
E = 0;

import casadi.*

% states
T = SX.sym('T', 1, 1); % torque (symbolic)
x = SX.sym('x',1,1);   % position
theta = SX.sym('theta',1,1); % wheel angle
u = SX.sym('u',1,1);         % velocity
omega = SX.sym('omega',1,1); % wheel angluar velocity
s = [x; theta; u; omega];    % states vector

% slip
% k =  -(u - omega.*r)./ if_else_smooth(u, 0, u+0.00001, -u+0.00001, 'C', 10);
k = -(u - omega.*r)./(abs(u) + 0.00001);

% M-F
y =  D.*sin( C.* atan( B.*k - E.*( B.*k - atan(B.*k))) );
y_fun = Function('y', {s}, {y});

% longitudinal force
Fx = if_else_smooth(   abs(T),...                                                % expression to check
                            0,...                                                % greater than
                    min([abs(y).*m.*9.81, abs(T./r), 1.2.*m.*9.81]).*sign(k),... % true
                    y.*m.*9.81,...                                               % false
                    'C', 20);                                                     % smooth parameter
Fx_fun = Function('Fx', {[s; T]}, {Fx});

% wheel dynamics 
% Euler equation
omegadot =  -(T + Fx.*r)./(0.5.*m.*r^2) ;
omegadot_fun = Function('od', {s, T}, {omegadot});

% Newton equation
Udot =  Fx./m;

% state-space dynamics
s_dot = [u; omega; Udot; omegadot];
s_dot_fun = Function('sd', {s, T}, {s_dot});

% initial condition
s_0 = [0;0;0.1;0]; % [x, theta, u, omega]


% simulation parameters
t0 = 0;
tf = 10;
dt_mean = 1e-3;
N = length(0:dt_mean:tf);
[~, ff] = RK4(s, T, s_dot, 5, tf, s_0, t0, zeros(1,N));

% INPUTS (torque values over the time horizon)
N = 1000;
nx = 4;
nu = 1;
nz = 1;
nuc = 0;
nzc = 0;
Npg = 1;
Xb.lb = [-inf; -inf; -inf; -inf];
Xb.ub = [inf; inf; inf; inf];
Xb.x0 = s_0;
Zb.lb = 1e-5;
Zb.ub = 1e-2;
Ub.lb = 0;
Ub.ub = inf;
Wgb.lb = 0;
Wgb.ub = 0;
d = 2;
Pdata.Pg = zeros(1,N); 
Pdata.Pg(1:floor(N/10)) = 0;
pb = nlp_casadi_2(nx, nu, nz, nuc, nzc, 0, N, d, Xb, Ub, Zb, [], [], [], Npg, 0, Pdata, 1, 0, 'SX');
pb.w0 = 0;

for i = 1:d
   pb.append_g(s_dot_fun(pb.xc(:,i), pb.u)*pb.z - pb.xp(:,i), zeros(pb.nx,1), zeros(pb.nx,1)) ;
end
% pb.append_g(s_dot_fun(pb.x, pb.u)*pb.z + pb.x_1 - pb.x, zeros(nx,1), zeros(nx,1));
% pb.Jj = (s_dot_fun(pb.x, pb.u)*pb.z)'*(s_dot_fun(pb.x, pb.u)*pb.z);
pb.append_g(pb.xc_end - pb.x , zeros(nx,1), zeros(nx,1));
pb.append_g(pb.u - pb.pg, 0, 0);

pb.append_g(sum(pb.Z) - tf, 0, 0, 'manually');
% Build map function 
pb.build_map;
% Build g 
pb.build_g;
% pb.build_J;

% integration
pb.build_solver([],'integrator')

% Solve NLP
sol_sym       = pb.solver('x0', pb.w0, 'lbx', pb.lbw, 'ubx', pb.ubw,'lbg', pb.lbg, 'ubg', pb.ubg);
siz           = length(sol_sym.x);  

x             = full(sol_sym.x(1:(nx + nu + nz +(nx)*d):siz))';
teta          = full(sol_sym.x(2:(nx + nu + nz +(nx)*d):siz))';
u_num         = full(sol_sym.x(3:(nx + nu + nz +(nx)*d):siz))';
omega         = full(sol_sym.x(4:(nx + nu + nz +(nx)*d):siz))';
u_in          = full(sol_sym.x(5:(nx + nu + nz +(nx)*d):siz))';
u_in(end+1)   = u_in(end);
h_opt         = full(sol_sym.x(6:(nx + nu + nz +(nx)*d):siz))';
t_num         = [0,cumsum(h_opt)];
s_sol         = [x; teta; u_num;omega];
%% animation
% figure; hold on; axis equal
% alpha = linspace(0, 2*pi, 50);
% pR = r.*[sin(alpha); cos(alpha)];
% tr = hgtransform(gca);
% plot(pR(1,:), pR(2,:), 'parent', tr);
% plot(pR(1,1), pR(2,1), 'marker', 'o', 'parent', tr)
% 
% x = full(s_sol(1, :));
% theta = full(s_sol(2, :));
% plot([x(1), x(end)], [0,0], 'color' , 'k', 'linewidth', 1.4); 
% for i = 1:10:length(x)
%    tr.Matrix = [rotZ(-theta(i)), [x(i); r; 0]; 0, 0, 0, 1];
%    xlim([x(i)-2, x(i) + 2])
%    pause(0.01)
% end


%% solution plots
% t_num = t0:dt:tf;
Fx_num = full(Fx_fun([s_sol; u_in]));
% omega = full(s_sol(4, :));
% u_num = full(s_sol(3, :));
k_num = -(u_num - omega.*r)./(abs(u_num)+0.0001);
slipping = (abs(u_num) - abs(omega.*r));
omega_dot_num = full(omegadot_fun(s_sol, u_in));

results = {Fx_num, omega, u_num, k_num, omega_dot_num, slipping};
labels = {'$F_x$','$\omega$','$u$','$k$','$\dot \omega$', 'slip'};
figure; hold on;
for i = 1:length(results)
    ax = subplot(3,2,i);
    plot(ax, t_num, results{i})
    xlabel('$t$ (s)');
    ylabel(labels{i});
    set(gca, 'FontSize', 22)
end
figure(4)
plot(t_num(1:end-1),h_opt)


