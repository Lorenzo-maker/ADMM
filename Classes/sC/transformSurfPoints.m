function [X, Y, Z] = transformSurfPoints(T, X, Y, Z)
% 
% apply rotation matrix to non vectorized set of points, such as surface
% grids
% [X, Y, Z] = rotateSurfPoints(R, X, Y, Z)
%

X = T(1,1).*X + T(1,2).*Y + T(1,3).*Z + T(1,4);
Y = T(2,1).*X + T(2,2).*Y + T(2,3).*Z + T(2,4);
Z = T(3,1).*X + T(3,2).*Y + T(3,3).*Z + T(3,4);
end