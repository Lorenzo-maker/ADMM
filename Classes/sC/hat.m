function X = hat(vec)
%
% trasforms a column vector R3x1 in antisym matrix form R3x3
%
%
% OR trasforms a twist column vector R6x1 in to the homogeneous
% hat form R4x4
%

if size(vec, 1) == 3
    X = [0 -vec(3) vec(2);...
        vec(3) 0 -vec(1);...
        -vec(2) vec(1) 0];
    return
end

X = [hat([vec(4);vec(5);vec(6)]), vec(1:3);0 0 0 0];


end