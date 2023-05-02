function Js = EulParSpatialJac(q)
%
% Js = EulParSpatialJac(q)
%
% returns the spatial jacobian corresponding to the Euler parameters such
% that :
% w^s = Js(q)*qdot
%
% input q: 4x1 array with the Euler parameters = [q0;q1;q2;q3]
%
q0 = q(1);     % scalar parameter
qvec = q(2:4); % vector part of the quaternion
Js = 2.*[-qvec, eye(3).*q0 + hat(qvec)];

end