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
nz = 1;
nuc = 0;
nzc = 0;
Xb.lb = [-inf; -inf];
Xb.ub = [inf; inf];
Xb.x0 = [1;1];
Zb.lb = 1e-5;
Zb.ub = 1e-2;
pb = nlp_casadi_2(nx,nu,nz,nuc,nzc,0,N,d, Xb, [], Zb, [], [], [], 0, 0, [], 1, 0, 'SX');

for i = 1:d
   pb.append_g(F(pb.xc(:,i))*pb.z - pb.xp(:,i), zeros(pb.nx,1), zeros(pb.nx,1)) ;
end
%
pb.append_g(pb.xc_end - pb.x , zeros(nx,1), zeros(nx,1));
pb.append_g(sum(pb.Z) - T, 0, 0, 'manually');
% Build map function 
pb.build_map;
% Build g 
pb.build_g;


% problem construction
prob = struct('f', pb.J, 'x', pb.w, 'g', pb.g);
opts = struct;
opts.ipopt.linear_solver = 'ma57';
opts.ipopt.max_iter = 5000;
opts.ipopt.print_level = 5;
solver = nlpsol('solver', 'ipopt', prob, opts);


% Solve NLP
sol_sym       = solver('x0', 0, 'lbx', pb.lbw, 'ubx', pb.ubw,'lbg', pb.lbg, 'ubg', pb.ubg);
siz           = length(sol_sym.x);  
sol.x_opt     = full(sol_sym.x(1:(nx + nz +(nx)*d):siz));
sol.y_opt     = full(sol_sym.x(2:(nx + nz +(nx)*d):siz));
sol.h_opt     = full(sol_sym.x(3:(nx + nz +(nx)*d):siz));
plot(sol.x_opt, sol.y_opt)



