function M = twistPole(d)
%
%
%
if length(d) > 3
    d = d(1:3,4);
end
cl = class(d);
I = eye(3, cl);

M = [I, hat(d);zeros(3,3, cl), I];

end