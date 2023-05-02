function [X, Y, Z] = rotateSurfPoints(R, X, Y, Z)
% 
% apply rotation matrix to non vectorized set of points, such as surface
% grids
% [X, Y, Z] = rotateSurfPoints(R, X, Y, Z)
%

X = R(1,1).*X + R(1,2).*Y + R(1,3).*Z;
Y = R(2,1).*X + R(2,2).*Y + R(2,3).*Z;
Z = R(3,1).*X + R(3,2).*Y + R(3,3).*Z;
end