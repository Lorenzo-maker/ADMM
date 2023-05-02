function J = floatingBaseSpatialJacQ(q)
%
% J = floatingBaseSpatialJacQ(q)
%
% returns the spatial jacobian of a floating base with the orientation
% described by Euler parameters
%
x = q(1:3); % positional variables
p = q(4:7); % orientation variables (Euler parameters)
G = EulParSpatialJac(p);
J = [eye(3), hat(x)*G; zeros(3,3, class(q(1))), G];
end