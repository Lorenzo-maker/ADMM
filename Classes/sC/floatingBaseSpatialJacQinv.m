function Jinv = floatingBaseSpatialJacQinv(q)
%
% J = floatingBaseSpatialJacQ(q)
%
% returns the spatial jacobian of a floating base with the orientation
% described by Euler parameters
%
x = q(1:3); % positional variables
p = q(4:7); % orientation variables (Euler parameters)
L = EulParSpatialJacInv(p);
Jinv = [eye(3), -hat(x); zeros(4,3, class(q(1))), L];
end