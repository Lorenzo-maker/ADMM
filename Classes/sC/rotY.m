function R = rotY(x)
%
c = cos(x);
s = sin(x);

R = [c 0 s;...
     0 1 0;...
    -s 0 c];
end