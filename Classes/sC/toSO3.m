function RinSO3 = toSO3(R)
% RinSO3 = toSO3(R)
% returns the best approximation of R in SO(3)
%

RinSO3 = nan(size(R));
S = eye(3);   % starting template for new S is always I
    
for i = 1:size(R, 3)
    [U, ~, V] = svd(R(:,:,i));
    S(3,3) = det(U*V');
    RinSO3(:,:,i) = U*S*V';
end

end
