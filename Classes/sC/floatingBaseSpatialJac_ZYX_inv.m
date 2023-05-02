function Jacinv = floatingBaseSpatialJac_ZYX_inv(q)

x = q(1:3);
alpha = q(4:6);
cl = class(q);
Jwinv = spatialJac_ZYX_inv(alpha);
Jacinv = [   eye(3, cl), -hat(x);...
          zeros(3,3,cl),  Jwinv];

end