function E = expTw(unitTwist,th)
%
% Exponential twist computation
%
% unitTwist = unitary twist of the joint
%
% th = joint variable value
%


if abs(unitTwist(4)) <= 1e-6 && abs(unitTwist(5)) <= 1e-6 && abs(unitTwist(6)) <= 1e-6 % prismatic joint
    
    axis = [unitTwist(1);unitTwist(2);unitTwist(3)];
    R = eye(3);
    d = axis*th;
    
else % revolute joint
    
    axis = [unitTwist(4);unitTwist(5);unitTwist(6)];
    v = [unitTwist(1);unitTwist(2);unitTwist(3)];
    q = -cross(v,axis);
    R = expSkew(axis,th);
    d =(eye(3) - R)*q;
    
end

E = [R d; 0 0 0 1];

end