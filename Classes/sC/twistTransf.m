function L = twistTransf(g)
%
% L = twistTransf(g)
%
% writes the matrix that maps the twist components to the new frame defined
% by the homogeneous transformation g or its rotation matrix
%
if size(g,2) == 4
    R = g(1:3, 1:3);
else
    R = g;
end
cl = class(g);
O = zeros(3, 3 , cl);
L = [R, O; O, R];


end