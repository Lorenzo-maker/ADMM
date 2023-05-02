function out = TrotZ(x)

c = cos(x);
s = sin(x);

out = [c,-s, 0, 0;...
       s, c, 0 ,0;...
       0, 0, 1, 0;...
       0, 0, 0, 1];

end