function out = TrotX(x)

c = cos(x);
s = sin(x);

out = [1, 0,  0, 0;...
       0, c, -s, 0;...
       0, s,  c, 0;...
       0, 0,  0, 1];

end