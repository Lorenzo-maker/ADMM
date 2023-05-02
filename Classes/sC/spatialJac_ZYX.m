function J = spatialJac_ZYX(alpha)

psi = alpha(1);
theta = alpha(2);

spsi = sin(psi);
cpsi = cos(psi);
ctheta = cos(theta);
stheta = sin(theta);

J = [0, -spsi, cpsi*ctheta;...
     0,  cpsi, ctheta*spsi;...
     1,     0,     -stheta];

end