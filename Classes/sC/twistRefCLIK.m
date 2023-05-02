function L = twistRefCLIK(taskM, eeM)
%
% reutrns matrix L to map the desired twist to the end-effector twist in
% order to obtain the dynamic orientation error vector parametrized with the
% geometric jacobian 
%
% inputs: taskM: homogeneous transform matrix of the task 
%           eeM: homogeneous transform matrix of end-effector

xd = taskM(1:3,1);
yd = taskM(1:3,2);
zd = taskM(1:3,3);
xee = eeM(1:3,1);
yee = eeM(1:3,2);
zee = eeM(1:3,3);
L = -0.5.*(hat(xd)*hat(xee) + hat(yd)*hat(yee) + hat(zd)*hat(zee));