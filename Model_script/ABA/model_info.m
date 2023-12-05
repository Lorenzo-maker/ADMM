disp('************************************************************************************')
disp('Model_info (remember to change info here and in car builder if something change)')
disp('States: [e_p, ef, d, theta, phi, alpha_dot, ep_dot, ef_dot, d_dot, theta_dot, phi_dot]')
nx = 11;
disp('Control: [Fz, Fx, Fy, Mx, My, Mz, delta]')
nu = 7;
disp('Algebraic: [Fzij, Yi, Fxa, Fxb]')
nz = 8;
disp('************************************************************************************')