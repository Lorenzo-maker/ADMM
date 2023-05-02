function out = rodriguezAxisAngle(ax, x)
%
% Rodriguez formula for the rotation matrix with axis angle
% parametrization
%

axisHat = hat(ax);
out = eye(3) + axisHat.*sin(x) + axisHat*axisHat.*(1 - cos(x));

end