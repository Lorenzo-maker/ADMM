function data = car_parameters_ocp_fun(alfarange, pista, lap, colloc_type, d)
% Computed grid
alfa_grid = alfarange;
pos_grid = full(pista.fun_pos(alfa_grid));
ns_grid  = full(pista.fun_vh(alfa_grid));
kp_grid  = full(pista.fun_kp(alfa_grid));
ts_grid = full(pista.fun_ts(alfa_grid));
ts_grid(3,:) = 0;
ts_grid = ts_grid./vecnorm(ts_grid);
psi_grid = atan_track(ts_grid, 'clockwise');
index0 = find(alfa_grid == 0);
if length(index0) > 1
    psi_laps = psi_grid(1:index0(2)-1);
    psi_lap = psi_grid(1:index0(2)-1);
    shift = psi_grid(index0(2)-1) - psi_grid(1);   
    for j = 2:length(index0)
        psi_laps = [psi_laps, psi_lap + (j-1)*shift];
    end
    psi_laps = [psi_laps, psi_lap(end) + (j-1)*shift];
    psi_grid = psi_laps;
end


load('Fsae_2022.mat')
[Fx0, Fy0, data.tyre.Fx, data.tyre.Fy, data.tyre.Fxlim, data.tyre.Fylim, data.tyre.Gxa, data.tyre.Gyk] = MF_Tyre(tir, 'shift', 0, 'curvature', 0);

% vehicle data
data.rho         = 1.2073;                    % air density
data.S           = 1.4;                       % front area
data.cx          = 0.84;                      % drag coefficient
data.cz          = 1.34;                      % lift coefficient
data.ab          = 0.4;
data.a1          = 0.765;
data.a2          = 0.815;
data.l           = data.a1 + data.a2;
data.hg          = 0.3;%0.28;
data.t1          = 1.2;
data.t2          = 1.2;
data.rw          = tir.UNLOADED_RADIUS;
data.m           = 240;                       % mass of car in kg
data.Jw          = (1+3)*0.5*tir.MASS*tir.UNLOADED_RADIUS^2;
data.Jz          = data.m*data.a1*data.a2*0.92;
data.G           = 9.81;                      % gravity
data.Fz_start    = data.m*data.G;             % weigth force (Newton)
data.mu_x        = 3.0 ;                      % dry road example
data.mu_y        = 3.0;                       % dry road example
data.tol         = 10^-3;                     % tolerance for avoiding NaN
data.Vm          = 22.23;                     % mean total speed (m/s) ~ 80 km/h
data.hm          = 0.3;                       % mean discretization 
data.Vi          = 30;                        % initial speed
data.Vmax        = 200/3.6;                   % max speed
data.Pmin        = -1800e3;                   % max braking power
data.Pmax        = 47e3;                      % max power
data.kb          = 0.6;
data.q1          = (2/5)*tir.UNLOADED_RADIUS;%0.06;
data.q2          = (1/5)*tir.UNLOADED_RADIUS;%0.065;
data.qm          = (data.a2*data.q1 + data.a1*data.q2)/data.l;
data.k_phi1      = 0.5*36000*data.t1^2;
data.k_phi2      = 0.5*24000*data.t2^2;
data.k_phi       = (data.k_phi1 + data.k_phi2 + 1e-5);
data.eta         = 1;
data.tau = 1/3;
data.d0 = 0;
data.ep1 = 0;

% Scaling Factors
data.u_scale  = data.Vmax;
data.v_scale  = 10;
data.r_scale  = 8;
data.x_max = max(pos_grid(1,:))+100;
data.x_min = min(pos_grid(1,:))-100;
data.x_scale = max(abs(pos_grid(1,:)));
data.y_max = max(pos_grid(2,:))+100;
data.y_min = min(pos_grid(2,:))-100;
data.y_scale = max(abs(pos_grid(2,:)));
data.psi_scale = 6*pi + max(lap-2,0)*2*pi;
data.w_scale = data.Vmax/data.rw;
data.Ta_scale = 2000;%data.Pmax/(data.w_scale);
data.Tb_scale = 10000;%data.Pmin/(data.w_scale);
data.delta_scale = pi;
data.Fz_scale = data.m*data.G + 0.5*data.rho*data.S*data.cz*data.Vmax^2;
data.Fx_scale = data.Fz_scale*data.mu_x;
data.Fy_scale = data.Fz_scale*data.mu_y;
data.P_scale  = max(abs(data.Pmin), data.Pmax);
data.ep_scale = 10;
data.h_scale = 1;
data.X_scale = [data.u_scale; data.v_scale; data.r_scale; data.x_scale; data.y_scale; data.psi_scale; repmat(data.w_scale,4,1)];
data.U_scale = [data.Ta_scale; data.Tb_scale; data.delta_scale];
data.Z_scale = [repmat(data.Fz_scale,4,1);repmat(data.Fx_scale,4,1);repmat(data.Fy_scale,4,1); data.ep_scale; data.h_scale];
X_scale = data.X_scale;
U_scale = data.U_scale;
Z_scale = data.Z_scale;
data.d = d;
data.colloc_type = colloc_type;
scale.x = X_scale;
scale.u = U_scale;
scale.z = Z_scale;
nx = length(data.X_scale);
nu = length(data.U_scale);
nz = length(data.Z_scale);
end


