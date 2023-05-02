clear all; clc; close all;
N = 100;
import casadi.*
% Funzioni
X = SX.sym('X',2,1);
U = SX.sym('U');
F = Function('F', {X, U}, {[X(2);U - X(2)]});
L = Function('L', {X(1)}, {1-sin(2*pi*X(1))/2});

Xb.lb = -inf*ones(2,1);
Xb.ub =  inf*ones(2,1);
Xb.x0 = [];
Ub.lb = -inf;
Ub.ub = inf;
Wgb.lb = -inf;
Wgb.ub = inf;
pb = nlp_casadi_2(2, 1, 0, 0, 0, 1, N, 0, Xb, Ub, [], [], [], Wgb, 0, 0, [], 1, 0, 'MX');

dt = pb.wg/N;
pb.append_g(pb.x_1 - [0;0], zeros(2,1), zeros(2,1), 'start')
% Integrator
k1 = F(pb.x_1,             pb.u);
k2 = F(pb.x_1 + dt/2*k1,   pb.u);
k3 = F(pb.x_1 + dt/2*k2,   pb.u);
k4 = F(pb.x_1 + dt*k3,     pb.u);
x_next = pb.x_1 + dt/6*(k1+2*k2+2*k3+k4); 
% Continuity
pb.append_g(pb.x - x_next, zeros(2,1), zeros(2,1)); % close the gaps
pb.append_g(pb.x(2) - L(pb.x(1)), -inf, 0);
pb.append_g(pb.u, 0, 1);
pb.append_g(pb.x(1) - 1, 0, 0, 'end');
pb.append_g(pb.wg, 0, inf, 'end');

pb.Jj = dt;
pb.cost;
pb.build_map;
pb.build_g;
pb.build_J;
pb.build_solver([], 'default')
pb.w0 = [repmat([0;1;0],N,1);[0;1;1]];
sol = pb.solver('x0', pb.w0, 'lbx', pb.lbw, 'ubx', pb.ubw,'lbg', pb.lbg, 'ubg', pb.ubg);
siz         = length(sol.x);  
sol.pos     = full(sol.x(1:3:siz));
sol.speed   = full(sol.x(2:3:siz));
sol.u       = full(sol.x(3:3:siz-1));

figure
hold on
plot(sol.speed);
plot(sol.pos);
plot(full(L(sol.pos)),'r--')
stairs(1:N,sol.u,'k');

sol.x(end)

