function out = smooth_max(x, options)
%
%
%
arguments
   x
   options.dim = 1;
end
dim = options.dim;
alpha = 5;
out = sum(x.*exp(alpha.*x), dim)./sum(exp(alpha.*x), dim);

end