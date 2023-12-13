function T = rtmInv(R,d)

T = [R.',-R.'*d;zeros(1,3),1];

end