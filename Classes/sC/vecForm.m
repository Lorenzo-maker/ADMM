function X = vecForm(Tw)
%
% transforms a hat form twist into the vec form twist: R_4x4 -> R_6
%
% or
%
% transforms a hat form antisymmetric matrix into the vector form: R_3x3 -> R_3

if size(Tw, 1) == 3
    
    X = [Tw(3,2);...
         Tw(1,3);...
         Tw(2,1)];
    return
end

w_hat = Tw(1:3,1:3);
v = Tw(1:3,4);
X = [v; vecForm(w_hat)];

end