function [n, theta] = rotToAxisAngle(R)

theta = acos((trace(R)-1)/2) ;
n = 1/(2*sin(theta)).*[R(3,2) - R(2,3); R(1,3) - R(3,1); R(2,1) - R(1,2)];

end