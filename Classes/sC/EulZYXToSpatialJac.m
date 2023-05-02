function out = EulZYXToSpatialJac(phi, theta, psi)
%
% returns the Spatial Jacobian Js corresponding to the ZYX Euler angles.
% This is such that w_s = Js dot{alpha, beta, gamma}, with w_s spatial \
% components
%
sphi = sin(phi);
cphi = cos(phi);
ctheta = cos(theta);
		out =  [0,    -sphi,            ctheta*cphi;...
                0,     cphi,            ctheta*sphi;...
                1,        0,            -sin(theta) ];
         
end