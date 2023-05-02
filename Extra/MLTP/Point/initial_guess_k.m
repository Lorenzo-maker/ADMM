import casadi.*

%% Initial Guess for States
u_0   = data.Vi;
v_0   = 0;
r_0   = 0;
x_0   = pos_grid(1,k+1);
y_0   = pos_grid(2,k+1);
psi_0 = psi_grid(k+1);

states_0        = [u_0, v_0, r_0, x_0, y_0, psi_0];
scaled_states_0 = (states_0./data.X_scale')';

%% Initial Guess for Controls

Fx_0 = (0.5*data.cx*data.S*data.rho*(u_0)^2);
rp_0 = 0;

controls_0        = [Fx_0, rp_0];
scaled_controls_0 = (controls_0./data.U_scale')';

%% Initial Guess for Algebraic Parameters

Fz_0 = data.m*data.G + 0.5*data.rho*data.S*data.cz*(u_0)^2;
Fy_0 = 0;
P_0  = Fx_0*u_0;
ep_0 = 0;
h_0  = 1;

algebraic_0        = [Fz_0, Fy_0, P_0, ep_0, h_0];
scaled_algebraic_0 = (algebraic_0./data.Z_scale')'; 
w0_k = [scaled_controls_0;scaled_algebraic_0; repmat(scaled_states_0, d, 1); scaled_states_0]; 

