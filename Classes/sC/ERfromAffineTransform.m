function [E, R] = ERfromAffineTransform(F)
%
% [E, R] = ERfromField(F) 
%
% extracts the Green deformation matrix (E) and the rigid rotation matrix (R)
% from the affine transform F
%
% It is just an alternative way of computing SVD for 3x3 matrices
%
B = F*F';
E = 0.5.*(B - eye(3));

[V, D] = eig(B);
D = sqrt(D);
invSqr = zeros(3);
for i = 1:3
    invSqr = invSqr + 1./D(i,i).*(V(:, i)*V(:, i)');
end

R = invSqr*F;

end