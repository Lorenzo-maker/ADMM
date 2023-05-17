%% Model Parameters (Build Model Structure)
%
model = fnc_vehicleParameterSet(pista, alpha_vec, Nsteps, d, lap, 0);
%
%% Car model with states, inputs and algebric parameters
[sys, rel] = fnc_vehicleModel(model, pista);

nx = length(sys.states.X);
nu = length(sys.inputs.U);
nz = length(sys.parameters.Z);

X_scale = sys.states.X_scale;
U_scale = [1;1;sys.inputs.delta_v_scaling];
Z_scale = [sys.parameters.Fxf_scaling; sys.parameters.Fxr_scaling; sys.parameters.Fyf_scaling;...
           sys.parameters.Fyr_scaling;sys.parameters.Fzf_scaling;sys.parameters.Fzr_scaling;...
           sys.parameters.ax_scaling;sys.parameters.ay_scaling;1;sys.parameters.gamma2_scaling;...
           sys.parameters.ep_scaling;sys.parameters.dt_scaling];
scale.x = X_scale;
scale.u = U_scale;
scale.z = Z_scale;       