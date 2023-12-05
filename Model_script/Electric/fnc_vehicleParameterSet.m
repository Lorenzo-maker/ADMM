function model = fnc_vehicleParameterSet(track, alfarange, N, d, lap, J_Homotopy, options)
    arguments
        track; alfarange; N; d; lap; J_Homotopy;
        options.alpha_in = 0; options.alpha_end = 1; 
    end
%
% Unit conversion factors
kph2mps = 3.6^-1;
kW2W = 1e3;
PS2W = 735.49875;
kWh2Ws = 3600 * kW2W;
C2K = 273.0;
rpm2radps = pi / 30;
%
% Global settings
model.global.g = 9.80665; % Gravity acceleration [m/s^2]
% Grid and Track evaluation
model.opt.d_colloc = d;
model.opt.dpts = N;
model.opt.lap = lap;
model.opt.J_Homotopy = J_Homotopy;
%model.track.cntrlpts = track.n;
model.track.alfa_grid = alfarange;
%alfa-grid
% if ~isempty(alfa_new)
%     model.track.alfa_grid_lap = alfa_new;
% else
%     model.track.alfa_grid_lap = linspace(options.alpha_in, options.alpha_end, model.opt.dpts/lap+1);
% end
% 
% if lap > 1
%     model.track.alfa_grid = [repmat(model.track.alfa_grid_lap(1:end-1),1,lap-1), model.track.alfa_grid_lap];
% else
%     model.track.alfa_grid = model.track.alfa_grid_lap;
% end
%track-quantities
model.track.pos_grid = full(track.fun_pos(model.track.alfa_grid));
model.track.ns_grid  = full(track.fun_vh(model.track.alfa_grid));
model.track.kp_grid  = full(track.fun_kp(model.track.alfa_grid));
model.track.ts_grid = full(track.fun_ts(model.track.alfa_grid));
model.track.ts_grid(3,:) = 0;
model.track.ts_grid = model.track.ts_grid./vecnorm(model.track.ts_grid);

%psi-track
if lap > 1
    model.track.psi_grid_lap = atan_track_diff(model.track.ts_grid(:,end - model.opt.dpts/lap  : end), 'clockwise');
    shift = model.track.psi_grid_lap(model.opt.dpts/lap) - model.track.psi_grid_lap(1);
    model.track.psi_grid = model.track.psi_grid_lap(1:model.opt.dpts/lap);
    for i = 1:lap-1
        model.track.psi_grid = [model.track.psi_grid, i*shift + model.track.psi_grid_lap(1:model.opt.dpts/lap)];       
    end
    model.track.psi_grid = [model.track.psi_grid, i*shift + model.track.psi_grid_lap(end)];
else
    model.track.psi_grid = atan_track_diff(model.track.ts_grid, 'clockwise');
end
model.track.RVSG = cell(1,model.opt.dpts+1);
for i = 1:model.opt.dpts+1
    model.track.RVSG{i} = [cos(-model.track.psi_grid(i)), sin(-model.track.psi_grid(i)), 0; -sin(-model.track.psi_grid(i)), cos(-model.track.psi_grid(i)), 0; 0, 0, 1]*(full(track.fun_Rgs(model.track.alfa_grid(i))))';
end
%
% Vehicle mass and inertia
model.vehicle.body.m = 2318.0; % Mass [kg]
model.vehicle.body.Jz = 3455.9; % Yaw inertia [kg-m^2]
%
% Vertical loads
model.vehicle.body.wbf = (1 - 1.4776/3.0021); % Weight balance front
model.vehicle.body.l = 3.0021; % Wheelbase [m]
model.vehicle.body.h = 0.5006; % Height of center of mass [m]
model.vehicle.body.Fxb_max = -2 * [9502; 4689]; % Maximum braking force of the axles [N]
model.vehicle.body.tau = [0.0693; 0.0]; % Axle steering reduction ratio
model.vehicle.body.C = model.vehicle.body.l * [(1 - model.vehicle.body.wbf); -model.vehicle.body.wbf]; % Axle center longitudinal coordinates relative to vehicle CG [m]
model.vehicle.body.Fz_0 = model.vehicle.body.m * model.global.g * [model.vehicle.body.wbf; (1 - model.vehicle.body.wbf)];
%
% Tire parameters
model.vehicle.tire.R_dyn = [0.354; 0.354]; % Tire dynamic radius [m]
model.vehicle.tire.mux_0 = [1.3706; 1.3706]; % Longitudinal coefficient of friction at zero load
model.vehicle.tire.k_mux = [0; 0]; % Load sensitivity of longitudinal coefficient of friction [N^-1]
model.vehicle.tire.muy_0 = [1.0933; 1.4027]; % Lateral coefficient of friction at zero load
model.vehicle.tire.k_muy = [0; 0]; % Load sensitivity of lateral coefficient of friction [N^-1]
model.vehicle.tire.By = [9.8252; 24.9526];
model.vehicle.tire.Cy = [1.6711; 0.7370];
model.vehicle.tire.Ey = [0.3190; -0.2693];
%
% Other parameters
model.track.width = 5; % Width of track [m]
%
% Aero parameters
model.vehicle.aero.rho = 1.225; % Air density [kg / m^3]
model.vehicle.aero.C_d = 0.317; % Drag coefficient
model.vehicle.aero.S = 2.656; % Frontal surface area [m^2]
model.vehicle.aero.C_l = [0.06; -0.05]; % Axle lift coefficient
model.vehicle.aero.h_d = 0.0; % Drag lever arm from ground [m]
%
model.opt.V_lower = 10;
model.opt.V_homotopy = 20;
% Performance parameters
model.opt.u_max = (kph2mps) * 400.0; %300.0; % Max longitudinal speed [m/s]
model.opt.v_max = 10.0; %5.0; % Maximum lateral speed [m/s]
model.opt.r_max = 6;
model.opt.x_max = max(abs(model.track.pos_grid(1,:)))+10; %maximum x-coordinate of the track centerline
model.opt.y_max = max(abs(model.track.pos_grid(2,:)))+10;
model.opt.psi_max = max(1,(lap/2))*6*pi;
model.opt.SOC_final = 0.15; % Final value of SOC greater than...
model.opt.T_max = (C2K) + 150.0; % Maximum temperature of the battery [K]
model.opt.ax_max = 20; % Maximum longitudinal acceleration [m/s^2]
model.opt.ay_max = 20; % Maximum lateral acceleration [m/s^2]
model.opt.ep_max = max(model.track.pos_grid(4,:))+0.1; % Maximum lateral distance [m]
dist = diff(model.track.pos_grid(1:2,:)');
dist = ceil(max(vecnorm(dist')));
model.opt.dt_max = 2*dist/(model.opt.V_lower);
% model.opt.dt_max = max(1,lap)*2000.0/model.opt.dpts; %1000.0; % Maximum time for the lap [s]
model.opt.delta_v_max = pi;
%
% other useful parameters
model.opt.V_avg = 30; % Mean total speed [m/s] [~ 80 km/h]
model.opt.u_0 = (kph2mps) * 175.0; %220.0; %100.0; % Initial speed [m/s]
model.opt.ep_0 = 0; % Initial lateral position [m]
model.opt.Fxt_0 = model.vehicle.aero.C_d * model.vehicle.aero.S ...
    * (0.5 * model.vehicle.aero.rho * model.opt.u_0^2); % Initial traction force [N]
model.opt.ep_margin = 0.85; % Minimum lateral distance to road borders [m]
%
% Implementation parameters
model.opt.gripxa_max = [1.00; 1.00]; % Maximum longitudinal grip utilization in acceleration
model.opt.gripxb_max = [1.00; 1.00]; % Maximum longitudinal grip utilization in braking
model.opt.gripy_max = [1.00; 0.90]; % Maximum lateral grip utilization
model.opt.tol = 1e-3; % Tolerance for avoiding NaN
model.opt.alpha_final = 0.99980; % Terminal condition on curve parameter
%
%%  Powertrain parameters
%
% e-Motor power curve parameters
model.vehicle.pt1.omega_base = (rpm2radps) * 4365.0;          
model.vehicle.pt1.P_max = (kW2W) * 160.0; % Power in constant sector [kW]
% e-Motor transmission
model.vehicle.pt1.eta = 0.96;
model.vehicle.pt1.tau_map_ratio = (10.65).^-1;
% model.vehicle.pt1.tau_map_speed = (kph2mps) * [999.9];
% model.vehicle.pt1.delta_speed = (kph2mps) * 10.0; % Speed difference required to reach 90% of next gear ratio
%
% e-Motor power curve parameters
model.vehicle.pt2.omega_base = (rpm2radps) * 6611.0;          
model.vehicle.pt2.P_max = (kW2W) * 270.0; % Power in constant sector [kW]
% e-Motor transmission
model.vehicle.pt2.eta = 0.96;
model.vehicle.pt2.tau_map_ratio = (10.65).^-1;
% model.vehicle.pt2.tau_map_speed = (kph2mps) * [999.9];
% model.vehicle.pt2.delta_speed = (kph2mps) * 10.0; % Speed difference required to reach 90% of next gear ratio
%
%%
% Battery parameters 
model.vehicle.bat.E_max = (kWh2Ws) * 84; % Nominal energy [Ws] 
model.vehicle.bat.eta_charge_tot = 0.87; % Charging efficiency with inverter
model.vehicle.bat.eta_discharge_tot = 0.87; % Discharging efficiency with inverter
model.vehicle.bat.Pcooling = (kW2W) * 11.0; % Power extracted by the cooling fluid [W]
model.vehicle.bat.COP = 3; % Coefficient of performance
model.vehicle.bat.eta_charge_batt = 0.90; % Charging efficiency only of the battery 
model.vehicle.bat.eta_discharge_batt = 0.90; % Discharging efficiency only of the battery
model.vehicle.bat.hc = 50; % Convective coefficient [W/(m^2-K)]
model.vehicle.bat.Ab = 4; % Area exposed of the battery [m^2]
model.vehicle.bat.mc = 533000; % Product of mass and specific heat coefficient [J/K]
model.vehicle.bat.T_0 = (C2K) + 25.0; % Initial battery temperature [K]
model.vehicle.bat.phi = 432 / (532 * 0.87 * 0.87); % Coefficient to consider the limit in charge for the battery
model.vehicle.bat.SOC_0 = 0.95; %0.95; % Initial SOC value
model.vehicle.bat.SOC_min = 0.10; % Min SOC value
model.vehicle.bat.SOC_max = 1.00; % Max SOC value
model.vehicle.bat.Pin_max = (kW2W) * 342.0; % Maximum Power in to the battery (regen) [kW]
%
end





