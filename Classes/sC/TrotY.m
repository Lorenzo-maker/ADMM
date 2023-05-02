function out = TrotY(x)
%
%
%
c = cos(x);
s = sin(x);


out = [c, 0, s, 0;...
       0, 1, 0, 0;...
      -s, 0, c, 0;...
       0, 0, 0, 1];

end