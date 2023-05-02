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
Xb.x0 = 'parametric';
Zb.lb = 1e-4;
Zb.ub = 1e-2;
Ub.lb = 0;
Ub.ub = inf;
Wgb.lb = 0;
Wgb.ub = 0;
np = nx + 1 + 1;
d = 2;
Pdata.Pg = zeros(1,N); 
Pdata.Pg(1:floor(N/10)) = 0;
pb = nlp_casadi_2(nx, nu, nz, nuc, nzc, 1, N, d, Xb, Ub, Zb, [], [], Wgb, Npg, 0, Pdata, 1, np, 'SX');
pb.w0 = 0;


for i = 1:d
   pb.append_g(s_dot_fun(pb.xc(:,i), pb.u)*pb.z - pb.xp(:,i), zeros(pb.nx,1), zeros(pb.nx,1)) ;
end
% pb.append_g(s_dot_fun(pb.x, pb.u)*pb.z + pb.x_1 - pb.x, zeros(nx,1), zeros(nx,1));
pb.append_g(pb.xc_end - pb.x , zeros(nx,1), zeros(nx,1));
% pb.Jj = (s_dot_fun(s_dot_fun(pb.x_1, pb.u)*(pb.z/2) + pb.x_1, pb.u)*(pb.z/2) + s_dot_fun(pb.x_1, pb.u)*(pb.z/2) + pb.x_1 -pb.x)'*(s_dot_fun(s_dot_fun(pb.x_1, pb.u)*(pb.z/2) + pb.x_1, pb.u)*(pb.z/2) + s_dot_fun(pb.x_1, pb.u)*(pb.z/2) + pb.x_1 -pb.x);
pb.append_g(pb.u - pb.p(5), 0, 0);

% Build map function 
pb.build_map;

% Add manually constraints that take in account the all variables (i.e. X, U, Z,...)
pb.append_g(sum(pb.Z) - pb.p(6), 0, 0, 'manually');

% Build g 
pb.build_g;

% pb.g = [pb.g; (sum(pb.Z) - pb.p(6))];%, 0, 0, 'manually');
% pb.lbg = [pb.lbg; 0];    
% pb.ubg = [pb.ubg; 0];    

pb.build_J;

% problem construction
pb.build_solver([], 'integrator');

% Solve NLP
N_steps = 10;
x = []; teta = []; u_num = []; omega = []; u_in = []; h_opt = []; opt_bool = [];
for i = 1:N_steps
    if i*(1/N_steps)*tf > 1
        torque = 0;
    else
        torque = 200;
    end
    if i == 1
        sol_sym       = pb.solver('x0', pb.w0, 'lbx', pb.lbw(s_0), 'ubx', pb.ubw(s_0),'lbg', pb.lbg, 'ubg', pb.ubg,'p', [s_0; torque; (1/N_steps)*tf]);
        siz           = length(sol_sym.x)-pb.nwg;  
        x           = [x, full(sol_sym.x(1:(nx + nu + nz +(nx)*d):siz))'];
        teta        = [teta, full(sol_sym.x(2:(nx + nu + nz +(nx)*d):siz))'];
        u_num       = [u_num, full(sol_sym.x(3:(nx + nu + nz +(nx)*d):siz))'];
        omega       = [omega, full(sol_sym.x(4:(nx + nu + nz +(nx)*d):siz))'];
        u_in        = [u_in, full(sol_sym.x(5:(nx + nu + nz +(nx)*d):siz))'];
        h_opt       = [h_opt, full(sol_sym.x(6:(nx + nu + nz +(nx)*d):siz))'];
        opt_bool = [opt_bool; pb.solver.stats.success];
    else
        s_1 = full(sol_sym.x(end- nx- pb.nwg+1:end - pb.nwg));
        sol_sym       = pb.solver('x0', sol_sym.x, 'lbx', pb.lbw(s_1), 'ubx', pb.ubw(s_1),'lbg', pb.lbg, 'ubg', pb.ubg,'p', [sol_sym.x(end-nx+1:end); torque; (1/N_steps)*tf]);
        siz           = length(sol_sym.x)-pb.nwg; 
        x           = [x, full(sol_sym.x((nx + nu + nz +(nx)*d) + 1:(nx + nu + nz +(nx)*d):siz))'];
        teta        = [teta, full(sol_sym.x((nx + nu + nz +(nx)*d) + 2:(nx + nu + nz +(nx)*d):siz))'];
        u_num       = [u_num, full(sol_sym.x((nx + nu + nz +(nx)*d) + 3:(nx + nu + nz +(nx)*d):siz))'];
        omega       = [omega, full(sol_sym.x((nx + nu + nz +(nx)*d) + 4:(nx + nu + nz +(nx)*d):siz))'];
        u_in        = [u_in, full(sol_sym.x(5:(nx + nu + nz +(nx)*d):siz))'];
        h_opt       = [h_opt, full(sol_sym.x(6:(nx + nu + nz +(nx)*d):siz))'];
        opt_bool = [opt_bool; pb.solver.stats.success];
    end
    
end
u_in(end+1) = u_in(end);
disp(opt_bool')
t_num       = [0,cumsum(h_opt)];
s_sol       = [x; teta; u_num;omega];
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
%    tr.Matrix = [rotZ(-teta(i)), [x(i); r; 0]; 0, 0, 0, 1];
%    xlim([x(i)-2, x(i) + 2])
%    pause(0.1)
% end


%% solution plots
Fx_num = full(Fx_fun([s_sol; u_in]));
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


