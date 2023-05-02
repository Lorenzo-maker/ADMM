function out = rotNthetaToQuat(n, theta)

p0 = cos(theta./2);
p1 = sin(theta./2).*n(1);
p2 = sin(theta./2).*n(2);
p3 = sin(theta./2).*n(3);

out = [p0;p1;p2;p3];
end