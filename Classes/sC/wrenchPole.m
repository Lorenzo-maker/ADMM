function W = wrenchPole(d)
%
%
%

I = eye(3);

W = [I, zeros(3,3, class(d));hat(d), I];

end