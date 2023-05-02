function R = rotX(alpha)

c = cos(alpha);
s = sin(alpha);

R = [1 0  0;...
     0 c -s;...
     0 s  c];
end