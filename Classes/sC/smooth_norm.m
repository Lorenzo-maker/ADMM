function out = smooth_norm(x)
%
%
%
if size(x, 2) > 1
    x  = x.';
end
out = sum(smooth_max([-x, x], 'dim', 2));
end