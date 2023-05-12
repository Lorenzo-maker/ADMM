
% Computed grid
alfa_grid = alpha_vec;
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

% if lap > 1
%     model.track.psi_grid_lap = atan_track(model.track.ts_grid(:,end - model.opt.dpts/lap  : end), track.type);
%     shift = model.track.psi_grid_lap(model.opt.dpts/lap) - model.track.psi_grid_lap(1);
%     model.track.psi_grid = model.track.psi_grid_lap(1:model.opt.dpts/lap);
%     for i = 1:lap-1
%         model.track.psi_grid = [model.track.psi_grid, i*shift + model.track.psi_grid_lap(1:model.opt.dpts/lap)];       
%     end
%     model.track.psi_grid = [model.track.psi_grid, i*shift + model.track.psi_grid_lap(end)];
% else
%     model.track.psi_grid = atan_track(model.track.ts_grid, track.type);
% end
%% Dati
data.rho         = 1.2073;                    % air density
data.S           = 1.4;                       % front area
data.cx          = 1.19;                      % drag coefficient
data.cz          = 3.63;                      % lift coefficient

data.G           = 9.81;                      % gravity
data.m           = 647.7;                     % mass of car in kg
data.Fz_start    = data.m*data.G;                       % weigth force (Newton)
data.mu_x        = 1.0 ;                      % dry road example
data.mu_y        = 1.0;                       % dry road example
data.tol         = 10^-3;                     % tolerance for avoiding NaN
data.Vm          = 22.23;                     % mean total speed (m/s) ~ 80 km/h
data.hm          = 0.25;                      % mean discretization 
data.Vi          = 10;                        % initial speed
data.Vmax        = 320/3.6;                   % max speed
data.Pmin        = -1800e3;                   % max braking power
data.Pmax        = 470e3;                     % max power

%% Scaling factors
data.u_scale  = data.Vmax;
data.v_scale  = 6;
data.r_scale  = 4;
data.x_max = max(pos_grid(1,:))+100;
data.x_min = min(pos_grid(1,:))-100;
data.x_scale = max(abs(pos_grid(1,:)));
data.y_max = max(pos_grid(2,:))+100;
data.y_min = min(pos_grid(2,:))-100;
data.y_scale = max(abs(pos_grid(2,:)));
data.psi_scale = 6*pi + max(lap-2,0)*2*pi;
data.Fz_scale = data.m*data.G + 0.5*data.rho*data.S*data.cz*data.Vmax^2;
data.Fx_scale = data.Fz_scale*data.mu_x;
data.rp_scale = 10;
data.Fy_scale = data.Fz_scale*data.mu_y;
data.P_scale  = max(abs(data.Pmin), data.Pmax);
data.ep_scale = 10;
data.h_scale = 5;
data.X_scale = [data.u_scale; data.v_scale; data.r_scale; data.x_scale; data.y_scale; data.psi_scale];
data.U_scale = [data.Fx_scale; data.rp_scale];
data.Z_scale = [data.Fz_scale; data.Fy_scale; data.P_scale; data.ep_scale; data.h_scale];
X_scale = data.X_scale;
U_scale = data.U_scale;
Z_scale = data.Z_scale;
scale.x = X_scale;
scale.u = U_scale;
scale.z = Z_scale;
nx = length(data.X_scale);
nu = length(data.U_scale);
nz = length(data.Z_scale);


