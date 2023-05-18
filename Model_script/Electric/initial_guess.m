%% Description of the file
% In this file we define our initial guesses for every nlp steprmset.
% This is made using a simple model (point mass model) and guessing the
% others quantities by it.

% Guessing position along the track
guess.u = model.opt.V_avg;
guess.v = 0;
%r guess
guess.r = kp_grid(k+1)*guess.u;
guess.x = pos_grid(1,k+1);
guess.y = pos_grid(2,k+1);
guess.psi = psi_grid(k+1);
% Guessing battery states
guess.SOC = model.vehicle.bat.SOC_0;
% DSOC = (model.opt.SOC_final - model.vehicle.bat.SOC_min) - model.vehicle.bat.SOC_0;
% guess.SOC = model.vehicle.bat.SOC_0 + k*DSOC/model.opt.dpts; 
guess.T = model.vehicle.bat.T_0;

% Guessing torque & brake split
guess.betaf = 1;
guess.betar = 1;
guess.gamma1 = 0;
guess.gamma2 = 0;

% Guessing input forces
guess.brk = 0;
guess.Fxt = model.vehicle.aero.C_d * model.vehicle.aero.S ...
    * (0.5 * model.vehicle.aero.rho * guess.u.^2);
guess.thr = 1;
% guess.thr = (guess.Fxt) * guess.u ...
%     / (rel.pt1_Pout_fun(guess.u / sys.states.u_scaling));

% Guessing vertical forces;
guess.Fzf =  model.vehicle.body.Fz_0(1) + model.vehicle.aero.C_l(1) ...
    * model.vehicle.aero.S * (0.5 * model.vehicle.aero.rho * guess.u^2);
guess.Fzr =  model.vehicle.body.Fz_0(2) + model.vehicle.aero.C_l(2) ...
    * model.vehicle.aero.S * (0.5 * model.vehicle.aero.rho * guess.u^2);

% Guessing longitudinal forces
guess.Fxf = guess.betaf ...
    .* rel.pt1_Pout_fun(guess.u / sys.states.u_scaling) ./ guess.u;
guess.Fxr = guess.betar ...
    .* rel.pt2_Pout_fun(guess.u / sys.states.u_scaling) ./ guess.u;

% Guessing accelerations
guess.ax = (guess.Fxf + guess.Fxr) ./ model.vehicle.body.m;
guess.ay = guess.u.^2 * kp_grid(k+1);

% Guessing lateral forces
guess.Fyf = (model.vehicle.body.Fz_0(1) / model.global.g) * guess.ay;
guess.Fyr = (model.vehicle.body.Fz_0(2) / model.global.g) * guess.ay;

% Guessing lateral displacement and dt
guess.ep = 0;
dist = diff(pos_grid(1:2,:)');
dist = ceil(mean(vecnorm(dist')));
guess.dt = dist/(guess.u);

% Guessing front slip angle
guess.alphayf = 0.05*if_else_smooth(guess.Fyf,0,1,-1,'C',0.00001);%sign(guess.Fyf);
% guess.alphayf = 0.05*sign(guess.Fyf);

% Guessing steer angle
guess.delta_v = (guess.alphayf + (guess.v ...
    + guess.r * model.vehicle.body.C(1)) ./ guess.u) ...
    / model.vehicle.body.tau(1);



% Switching to normalized values for states
guess.u_norm = (sys.states.u_scaling^-1) * guess.u;
guess.v_norm = (sys.states.v_scaling^-1) * guess.v;
guess.r_norm = (sys.states.r_scaling^-1) * guess.r;
guess.x_norm = (sys.states.x_scaling^-1) * guess.x;
guess.y_norm = (sys.states.y_scaling^-1) * guess.y;
guess.psi_norm = (sys.states.psi_scaling^-1) * guess.psi;
guess.SOC_norm = (sys.states.SOC_scaling^-1) * guess.SOC; 
guess.T_norm = (sys.states.T_scaling^-1) * guess.T;

% Switching to normalized values for inputs
guess.delta_v_norm = (sys.inputs.delta_v_scaling^-1) * guess.delta_v;

% Switching to normalized values for parameters
guess.Fxf_norm = (sys.parameters.Fxf_scaling^-1) * guess.Fxf;
guess.Fxr_norm = (sys.parameters.Fxr_scaling^-1) * guess.Fxr;
guess.Fyf_norm = (sys.parameters.Fyf_scaling^-1) * guess.Fyf;
guess.Fyr_norm = (sys.parameters.Fyr_scaling^-1) * guess.Fyr;
guess.Fzf_norm = (sys.parameters.Fzf_scaling^-1) * guess.Fzf;
guess.Fzr_norm = (sys.parameters.Fzr_scaling^-1) * guess.Fzr;
guess.ax_norm = (sys.parameters.ax_scaling^-1) * guess.ax;
guess.ay_norm = (sys.parameters.ay_scaling^-1) * guess.ay;
guess.gamma1_norm = (sys.parameters.gamma1_scaling^-1) * guess.gamma1;
guess.gamma2_norm = (sys.parameters.gamma2_scaling^-1) * guess.gamma2;
guess.ep_norm = (sys.parameters.ep_scaling^-1) * guess.ep;
guess.dt_norm = (sys.parameters.dt_scaling^-1) * guess.dt;

% Input vector guessing at step k
sys.inputs.U_0 = [ ...
    guess.brk; ...
    guess.thr; ...
    guess.delta_v_norm; ...
    ];

% State vector guessing at step k
sys.states.X_0 = [ ...
    guess.u_norm; ...
    guess.v_norm; ...
    guess.r_norm; ...
    guess.x_norm; ...
    guess.y_norm; ...
    guess.psi_norm; ...
    guess.SOC_norm; ...
    guess.T_norm; ...
    ];

% Parameters vector guessing at step k
sys.parameters.Z_0 = [ ...
    guess.Fxf_norm; ...
    guess.Fxr_norm; ...
    guess.Fyf_norm; ...
    guess.Fyr_norm; ...
    guess.Fzf_norm; ...
    guess.Fzr_norm; ...
    guess.ax_norm; ...
    guess.ay_norm; ...
    guess.betaf; ...
    guess.betar; ...
    guess.gamma1_norm; ...
    guess.gamma2_norm; ...
    guess.ep_norm; ...
    guess.dt_norm; ...
    ];
%
