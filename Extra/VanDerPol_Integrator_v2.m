clear all
clc
import casadi.*
X = SX.sym('X', 2, 1);
Xdot = [(1-X(2)^2)*X(1) - X(2); X(1)];
F = Function('F', {X}, {Xdot});
T = 10;
d = 2;
N = 1000;
nx = 2;
nu = 0;
nz = 0;
nuc = 0;
nzc = 0;
Xb.lb = [-inf; -inf];
Xb.ub = [inf; inf];
Xb.x0 = 'parametric';
pb = nlp_casadi_2(nx,nu,nz,nuc,nzc,0,1,d, Xb, [], [], [], [], [], 0, 0, [], 1, nx, 'SX');

% pb.append_g(pb.x_1 - pb.p , zeros(nx,1), zeros(nx,1), 'start');

for i = 1:d
   pb.append_g(F(pb.xc(:,i))*(T/N) - pb.xp(:,i), zeros(pb.nx,1), zeros(pb.nx,1)) ;
end
%
pb.append_g(pb.xc_end - pb.x, zeros(nx,1), zeros(nx,1));
% Build map function 
pb.build_map;
% Build g 
pb.build_g;

% problem construction
prob = struct('f', pb.J, 'x', pb.w, 'g', pb.g, 'p', pb.p);
opts = struct;
opts.ipopt.linear_solver = 'ma57';
opts.ipopt.max_iter = 5000;
opts.ipopt.print_level = 0;
opts.print_time = 0;
solver = nlpsol('solver', 'ipopt', prob, opts);

xx = [];
yy = [];
opt_bool = [];
% Solve NLP
for i = 1:N
    if i == 1
        s_0 = [1;1];
        sol_sym       = solver('x0', 0, 'lbx', pb.lbw(s_0), 'ubx', pb.ubw(s_0),'lbg', pb.lbg, 'ubg', pb.ubg, 'p', s_0);
        siz           = length(sol_sym.x);  
        xx = [xx,sol_sym.x(1)];
        yy = [yy,sol_sym.x(2)];
        opt_bool = [opt_bool; solver.stats.success];
    else 
        s_0 = sol_sym.x(end-1:end);
        sol_sym       = solver('x0', sol_sym.x, 'lbx', pb.lbw(s_0), 'ubx', pb.ubw(s_0),'lbg', pb.lbg, 'ubg', pb.ubg, 'p', s_0);
        siz           = length(sol_sym.x);  
        xx = [xx,sol_sym.x(1)];
        yy = [yy,sol_sym.x(2)];
        opt_bool = [opt_bool; solver.stats.success];        
    end
end
if (all(opt_bool == 1))
   disp('all problems converged')
end
plot(full(xx), full(yy))




