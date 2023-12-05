import casadi.*

%% Initial Guess for States
u_0   = car.data.Vi;
v_0   = 0;
r_0   = 0*kp_grid(k+1)*u_0;
x_0   = pos_grid(1,k+1);
y_0   = pos_grid(2,k+1);
psi_0 = psi_grid(k+1);
omega_0 = u_0/car.data.rw;


states_0        = [u_0, v_0, r_0, x_0, y_0, psi_0, repmat(omega_0, 1, 2), repmat(omega_0, 1, 2)];
scaled_states_0 = (states_0./car.data.X_scale')';

%% Initial Guess for Controls

Ta_0 = (0.5*car.data.cx*car.data.S*car.data.rho*(u_0)^2)*car.data.rw;
Tb_0 = 0;
delta_0 = 0*car.data.m*u_0*r_0/30e3;

controls_0        = [Ta_0, Tb_0, delta_0];
scaled_controls_0 = (controls_0./car.data.U_scale')';

%% Initial Guess for Algebraic Parameters

Fz_0 = 0.25*(car.data.m*car.data.G + 0.5*car.data.rho*car.data.S*car.data.cz*(u_0)^2);
Fx1_0 = 0;
Fx2_0 = Ta_0/2/car.data.rw;
Fy_0 = car.data.m*u_0*r_0/4;
ep_0 = 0;
h_0  = 0.01;

algebraic_0        = [repmat(Fz_0,1,4),Fx1_0,Fx1_0,Fx2_0,Fx2_0,repmat(Fy_0,1,4), ep_0, h_0];
scaled_algebraic_0 = (algebraic_0./car.data.Z_scale')'; 
w0_k = [scaled_controls_0;scaled_algebraic_0; repmat(scaled_states_0, d, 1); scaled_states_0]; 

