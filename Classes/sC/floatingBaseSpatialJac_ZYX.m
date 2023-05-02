function Jac = floatingBaseSpatialJac_ZYX(q)

x = q(1:3);
alpha = q(4:6);
cl = class(q);
Jw = spatialJac_ZYX(alpha);
Jac = [   eye(3, cl), hat(x)*Jw;...
       zeros(3,3,cl),        Jw];

end