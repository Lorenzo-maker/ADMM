function out = Euler_FS(x, u, xdot, dt, options)
%
% x_k+1 = Euler_FS(x, u, xdot, dt, options)
% Euler forward step of a dynamic system expressed in state space:
% xdot = f(x, u);
% x = state space variable vector
% u = input variable vector
% dt is the integration time integral
% 
% options:
%
% 'Nsteps' (default = 1) sets the number of integration steps per
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
    
    out = xdot(out, u).*h + out;
    
end

end
