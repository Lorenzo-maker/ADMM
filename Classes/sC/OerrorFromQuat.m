function e_o = OerrorFromQuat(Qtask, QEE)
%
%
%
p0T = Qtask(1);
pvecT = Qtask(2:4);

p0E = QEE(1);
pvecE = QEE(2:4);

e_o = -p0T.*pvecE + p0E.*pvecT - hat(pvecT)*pvecE;
end