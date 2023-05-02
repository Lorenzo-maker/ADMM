function [surfObj, X, Y, Z] = transformSurface(surfObj, T)
% 
% apply homogeneous transformation matrix to a surface
% [surfObj, XData, YData, ZData] = transformSurface(surfObj, T)
%
% returns surface object and the grid data

R = T(1:3, 1:3);
d = T(1:3, 4);
XData = surfObj.XData;
YData = surfObj.YData;
ZData = surfObj.ZData;
X = R(1,1).*XData + R(1,2).*YData + R(1,3).*ZData + d(1);
Y = R(2,1).*XData + R(2,2).*YData + R(2,3).*ZData + d(2);
Z = R(3,1).*XData + R(3,2).*YData + R(3,3).*ZData + d(3);

surfObj.XData = X;
surfObj.YData = Y;
surfObj.ZData = Z;


end