function max = soft_max(x)
%
%
%

max = log(sum(exp(10.*x)))./10;

end