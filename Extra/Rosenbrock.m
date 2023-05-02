clear all
close all
clc

import casadi.*
x = SX.sym('x');
y = SX.sym('y');
z = SX.sym('z');
X = [x;y;z];

f = x^2 + 100*z^2;
L = Function('L', {X}, {f});
Zb.lb = [-inf; -inf; -inf];
Zb.ub = [inf; inf; inf];
N = 1;
d = 0;
pb = nlp_casadi_2(0, 0, 3, 0, 0, 0, N, d, [], [], Zb, [], [], [], 0, 0, [], 1, 0, 'SX');

pb.Jj = L(pb.z);
pb.append_g(pb.z(3) - pb.z(2) + (1-pb.z(1))^2 , 0, 0);

pb.cost;
pb.build_map;
pb.build_J;
pb.build_g;

% problem construction
prob = struct('f', pb.J, 'x', pb.w, 'g', pb.g);
opts = struct;
opts.ipopt.linear_solver = 'ma57';
opts.ipopt.max_iter = 5000;
opts.ipopt.print_level = 5;
solver = nlpsol('solver', 'ipopt', prob, opts);


% Solve NLP
sol_sym       = solver('x0', repmat([2.5; 3; 0.75], N, 1), 'lbx', pb.lbw, 'ubx', pb.ubw,'lbg', pb.lbg, 'ubg', pb.ubg);
siz           = length(sol_sym.x);  
sol_sym.x
x_opt = sol_sym.x(1:(2+(2)*pb.d):siz);
y_opt = sol_sym.x(2:(2+(2)*pb.d):siz);