function R = rotZ(x)
%
%
%

c = cos(x);
s = sin(x);

R = [c -s  0;...
     s  c  0;...
     0  0  1];

end