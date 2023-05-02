function J = spatialJac_ZYX_inv(alpha)

psi = alpha(1);
theta = alpha(2);

spsi = sin(psi);
cpsi = cos(psi);
ctheta = cos(theta);
stheta = sin(theta);

J = [(cpsi*stheta)/ctheta, (spsi*stheta)/ctheta, 1;...
                    -spsi,                 cpsi, 0;...
              cpsi/ctheta,          spsi/ctheta, 0];

end