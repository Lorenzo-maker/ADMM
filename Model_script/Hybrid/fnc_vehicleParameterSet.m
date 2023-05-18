function model = fnc_vehicleParameterSet(track, alfa_new, N, d, lap, J_Homotopy, options)
    arguments
        track; alfa_new; N; d; lap; J_Homotopy;
        options.alpha_in = 0; options.alpha_end = 1; 
    end
%
%% JSN P4 Model Parameters (22.03.14)
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
%
% Grid and Track evaluation
model.opt.d_colloc = d;
model.opt.dpts = N;
model.opt.J_Homotopy = J_Homotopy;
model.track.cntrlpts = track.n;
%alfa-grid
if ~isempty(alfa_new)
    model.track.alfa_grid_lap = alfa_new;
else
    model.track.alfa_grid_lap = linspace(options.alpha_in,options.alpha_end,model.opt.dpts/lap+1);
end

if lap > 1
    model.track.alfa_grid = [repmat(model.track.alfa_grid_lap(1:end-1),1,lap-1), model.track.alfa_grid_lap];
else
    model.track.alfa_grid = model.track.alfa_grid_lap;
end
%track-quantities
model.track.pos_grid = full(track.fun_pos(model.track.alfa_grid));
model.track.ns_grid  = full(track.fun_vh(model.track.alfa_grid));
model.track.kp_grid  = full(track.fun_kp(model.track.alfa_grid));
model.track.ts_grid = full(track.fun_ts(model.track.alfa_grid));
model.track.ts_grid(3,:) = 0;
model.track.ts_grid = model.track.ts_grid./vecnorm(model.track.ts_grid);
%psi-track
if lap > 1
    model.track.psi_grid_lap = atan_track(model.track.ts_grid(:,end - model.opt.dpts/lap  : end), 'clockwise');
    shift = model.track.psi_grid_lap(model.opt.dpts/lap) - model.track.psi_grid_lap(1);
    model.track.psi_grid = model.track.psi_grid_lap(1:model.opt.dpts/lap);
    for i = 1:lap-1
        model.track.psi_grid = [model.track.psi_grid, i*shift + model.track.psi_grid_lap(1:model.opt.dpts/lap)];
    end
    model.track.psi_grid = [model.track.psi_grid, i*shift + model.track.psi_grid_lap(end)];
else
    model.track.psi_grid = atan_track(model.track.ts_grid, 'clockwise');
end
model.track.RVSG = cell(1,model.opt.dpts+1);
for i = 1:model.opt.dpts+1
    model.track.RVSG{i} = [cos(-model.track.psi_grid(i)), sin(-model.track.psi_grid(i)), 0; -sin(-model.track.psi_grid(i)), cos(-model.track.psi_grid(i)), 0; 0, 0, 1]*(full(track.fun_Rgs(model.track.alfa_grid(i))))';
end
% Vehicle mass and inertia
model.vehicle.body.m = 1779.2; % Mass [kg]
model.vehicle.body.Jz = 2998.9; % Yaw inertia [kg-m^2]
%
% Vertical loads
model.vehicle.body.wbf = 0.5574; % Weight balance front
model.vehicle.body.l = 2.650; % Wheelbase [m]
model.vehicle.body.h = 0.481; % Height of center of mass [m]
model.vehicle.body.Fxb_max = -2 * [7035; 3106]; % Maximum braking force of the axles [N]
model.vehicle.body.tau = [0.0721; 0.0]; % Axle steering reduction ratio
model.vehicle.body.C = model.vehicle.body.l * [(1 - model.vehicle.body.wbf); -model.vehicle.body.wbf]; % Axle center longitudinal coordinates relative to vehicle CG [m]
model.vehicle.body.Fz_0 = model.vehicle.body.m * model.global.g * [model.vehicle.body.wbf; (1 - model.vehicle.body.wbf)];
%
% Tire parameters
model.vehicle.tire.R_dyn = [0.312; 0.312]; % Tire dynamic radius [m]
model.vehicle.tire.mux_0 = [1.1262; 1.2894]; % Longitudinal coefficient of friction at zero load
model.vehicle.tire.k_mux = [0; 0]; % Load sensitivity of longitudinal coefficient of friction [N^-1]
model.vehicle.tire.muy_0 = [1.1262; 1.2894]; % Lateral coefficient of friction at zero load
model.vehicle.tire.k_muy = [0; 0]; % Load sensitivity of lateral coefficient of friction [N^-1]
model.vehicle.tire.By = [9.6013; 20.1161];
model.vehicle.tire.Cy = [1.5302; 0.9190];
model.vehicle.tire.Ey = [-0.3653; -0.1902];
%
% Other parameters
model.track.width = 5; % Width of track [m]
%
% Aero parameters
model.vehicle.aero.rho = 1.225; % Air density [kg / m^3]
model.vehicle.aero.C_d = 0.334; % Drag coefficient
model.vehicle.aero.S = 2.230; % Frontal surface area [m^2]
model.vehicle.aero.C_l = [0.00; 0.00]; % Axle lift coefficient
model.vehicle.aero.h_d = 0.0; % Drag lever arm from ground [m]
%
model.opt.V_lower = 10;
model.opt.V_homotopy = 20;
% Performance parameters
model.opt.u_max = (kph2mps) * 400.0; % Max longitudinal speed [m/s]
model.opt.v_max = 10.0; %5.0; % Maximum lateral speed [m/s]
model.opt.r_max = 6; % Maximum yaw rate speed
model.opt.x_max = max(abs(model.track.pos_grid(1,:)))+10; %maximum x-coordinate of the track centerline
model.opt.y_max = max(abs(model.track.pos_grid(2,:)))+10;
model.opt.psi_max = max(1,(lap/2))*6*pi;
model.opt.SOC_final = 0.55; % Final value of SOC greater than...
model.opt.T_max = (C2K) + 35; %150.0; % Maximum temperature of the battery [K]
model.opt.ax_max = 20; % Maximum longitudinal acceleration [m/s^2]
model.opt.ay_max = 20; % Maximum lateral acceleration [m/s^2]
model.opt.ep_max = max(model.track.pos_grid(4,:))+0.1; % Maximum lateral distance [m]
model.opt.time_max = max(1,lap/2)*2000.0; % Maximum time for the lap [s]
dist = diff(model.track.pos_grid(1:2,:)');
dist = ceil(max(vecnorm(dist')));
model.opt.dt_max = 2*dist/(model.opt.V_lower);
% model.opt.dt_max = max(1,lap)*2000.0/model.opt.dpts; %1000.0; % Maximum time for the lap [s]
model.opt.delta_v_max = pi;
model.opt.ddelta_v_max = pi;
%
model.opt.V_avg = 30; %16; % Mean total speed [m/s] [~ 80 km/h]
model.opt.u_0 = (kph2mps) * 220.0; %175.0; % Initial speed [m/s]
model.opt.ep_0 = 0; % Initial lateral position [m]
model.opt.Fxt_0 = model.vehicle.aero.C_d * model.vehicle.aero.S ...
    * (0.5 * model.vehicle.aero.rho * model.opt.u_0^2); % Initial traction force [N]
model.opt.ep_margin = 0.85;%2.85; %0.85; % Minimum lateral distance to road borders [m]

% Implementation parameters
model.opt.gripxa_max = [1.00; 1.00]; % Maximum longitudinal grip utilization in acceleration
model.opt.gripxb_max = [1.00; 1.00]; % Maximum longitudinal grip utilization in braking
model.opt.gripy_max = [1.00; 1.00]; % Maximum lateral grip utilization
model.opt.tol = 1e-3; % Tolerance for avoiding NaN
model.opt.alpha_final = 0.99980; % Terminal condition on curve parameter
%
%%  Powertrain parameters
%
% IC-Engine power curve parameters
model.vehicle.pt1.Pout_map_speed = (rpm2radps) * [1008; 1500; 2000; 2503; 3002; 3504; 4005; 4503; 5006; 5502; 6003; 6500; 6800];
model.vehicle.pt1.Pout_map_power = (PS2W) * [31.2; 73.9; 102.2; 127.6; 152.1; 175.3; 198.1; 224.2; 244.7; 259.3; 268.4; 263.2; 246.9];
model.vehicle.pt1.P_map = polyfit(model.vehicle.pt1.Pout_map_speed, model.vehicle.pt1.Pout_map_power, 2);
% IC-Engine transmission
model.vehicle.pt1.eta = 0.86; %0.735;
model.vehicle.pt1.tau_map_ratio = ([3.714; 2.261; 2.174; 1.621; 0.927; 0.767; 0.878; 0.698] .* ...
    [3.800; 3.800; 2.714; 2.714; 3.800; 3.800; 2.714; 2.714]).^-1;
model.vehicle.pt1.tau_map_speed = (kph2mps) * [52.6; 86.8; 126.4; 169.5; 211.8; 324.9; 397.4; 999.9];
model.vehicle.pt1.delta_speed = (kph2mps) * 10.0; % Speed difference required to reach 90% of next gear ratio
%
% e-Motor power curve parameters
model.vehicle.pt2.omega_base = (rpm2radps) * 4000.0;          
model.vehicle.pt2.P_max = (kW2W) * 83.0; % Power in constant sector [kW]
% e-Motor transmission
model.vehicle.pt2.eta = 1.00;
model.vehicle.pt2.tau_map_ratio = ([14.80; 5.05]).^-1;
model.vehicle.pt2.tau_map_speed = (kph2mps) * [102.0; 999.9];
model.vehicle.pt2.delta_speed = (kph2mps) * 10.0; % Speed difference required to reach 90% of next gear ratio
%
%%
% Battery parameters 
model.vehicle.bat.E_max = (kWh2Ws) * 11.1; % Nominal energy [Ws] 
model.vehicle.bat.eta_charge_tot = 0.857; % Charging efficiency with inverter
model.vehicle.bat.eta_discharge_tot = 0.857; % Discharging efficiency with inverter
model.vehicle.bat.Pcooling = (kW2W) * 3.0; % Power extracted by the cooling fluid [W]
model.vehicle.bat.COP = 3; % Coefficient of performance
model.vehicle.bat.eta_charge_batt = 0.9; % Charging efficiency only of the battery 
model.vehicle.bat.eta_discharge_batt = 0.9; % Discharging efficiency only of the battery
model.vehicle.bat.hc = 50; % Convective coefficient [W/(m^2-K)]
model.vehicle.bat.Ab = 4; % Area exposed of the battery [m^2]
model.vehicle.bat.mc = 50700; % Product of mass and specific heat coefficient [J/K]
model.vehicle.bat.T_0 = (C2K) + 30; %25.0; % Initial battery temperature [K]
model.vehicle.bat.SOC_0 = 0.75; % Initial SOC value
model.vehicle.bat.SOC_min = 0.10; % Min SOC value
model.vehicle.bat.SOC_max = 0.95; % Max SOC value
model.vehicle.bat.Pin_max = (kW2W) * 66.0; % Maximum Power in to the battery (regen) [kW]
%
end





