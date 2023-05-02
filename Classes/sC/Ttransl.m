function out = Ttransl(d)

out = eye(4, class(d));
out(1:3,4) = d;

end