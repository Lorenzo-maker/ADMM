function GstInv = rigidInverse(Gst)
%
% Computes the inverse of a homogeneous transormation Gst into
% Gst^-1


Rts = Gst(1:3,1:3).';
dst = Gst(1:3,4);
GstInv = [Rts,-Rts*dst; 0 0 0 1];
end