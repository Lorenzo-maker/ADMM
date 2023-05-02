function [xNext, S] = Euler_BS(x, u, xdot, dt, options)
%
% x_k+1 = Euler_BS(x, u, xdot, dt, options)
% Euler backwards step of a dynamic system expressed in state space:
% xdot = f(x, u);
% x = state space variable vector (numeric)
% u = input variable vector       (numeric)
%
% Euler backwards step is an implicit method that solves the (in general)
% non linear equation:
% x_k+1 = x_k + h*f(x_k+1, t_k+1, u);
%
% dt is the integration time integral
% 
% options:
%  'Nsteps' (default = 1) sets the number of integration steps per
%           integration interval
%  'Solver' (default []) sets the nonlinear rootfining solver for the
%           implicit integration step. If the field is empty a new solver
%           is built

arguments
   x
   u
   xdot
   dt
   options.Nsteps = 1;
   options.Solver = [];
end

S = options.Solver;
Nsteps = options.Nsteps;

h = dt./Nsteps;
xNext = x;

if isempty(S) % build a casadi findroot solver once
    import casadi.*
    x2 = SX.sym('x2', length(x), 1); % x_k+1
    x1 = SX.sym('x1', length(x), 1); % x_k
    eq = x2 - xdot(x2, u).*h - x1;
    problem = struct('x', x2, 'p', x1, 'g', eq);
    S = rootfinder('S', 'newton', problem);

end

for i = 1:Nsteps
    
      xNext = S('x0', xNext, 'p', xNext);
      xNext = full(xNext.x);
end

end