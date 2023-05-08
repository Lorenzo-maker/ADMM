function out = RK4_step(x, u, xdot, dt, options)
%
% x_k+1 = RK4_step(x, u, xdot, dt, options)
% Euler forward step of a dynamic system expressed in state space:
% xdot = f(x, u) as a casadi Function;
% x = state space variable vector
% u = input variable vector
% dt is the integration time integral
% 
% options:
%  'Nsteps' (default = 1) sets the number of integration steps per
%           integration interval

arguments
   x
   u
   xdot
   dt
   options.Nsteps = 1;
end
Nsteps = options.Nsteps;

h = dt./Nsteps;
out = x;

% compute the integration through the time interval dt
for i = 1:Nsteps
    
    a1 = xdot(out           , u);
    a2 = xdot(out + h./2.*a1, u);
    a3 = xdot(out + h./2.*a2, u);
    a4 = xdot(out + h.*a3   , u);
    out = out + h./6.*(a1 + 2.*a2 + 2.*a3 + a4);
    
end

end