function out = if_else_smooth(in1, in2, T, F, options)
%
%
%
arguments
   in1
   in2
   T
   F
   options.C = 10;
end
gt = tanh(options.C.*(in1 - in2));
out = T./2.*(1 + gt) + (F./2).*(1-gt);
end