function R = quatToRot(Quat)
%
% R = quatToRot(Quat)
%
% returns the rotation matrix from the Euler parameters.
% Quat: 4x1 vector with all the Euler parameters = [q0;q1;q2;q3]

q0 = Quat(1);  % scalar parameter
q = Quat(2:4); % vector parameters

hq = hat(q);
I = eye(3);
R = I + 2*hq*(q0.*I + hq);


end