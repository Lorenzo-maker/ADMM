% Script for launching ad hoc car_builder
car = car_builder(alpha_vec, pista, colloc_type, d, lap, Nsteps, 0);

function car = car_builder(alfarange, pista, colloc_type, d, lap, Nsteps, sub_homotopy)
    import casadi.*
    %% Model Parameters (Build Model Structure)
    %
    model = fnc_vehicleParameterSet(pista, alfarange, Nsteps, d, lap, 0);
    %
    model.opt.colloc_type = colloc_type;
    %% Car model with states, inputs and algebric parameters
    [sys, rel] = fnc_vehicleModel(model, pista);
    
    driver.K1 = 0.95;%0.95;
    driver.K2 = 0.025;%0.025;
    driver.K3 = 0.025;%0.025;
    driver.der_r_factor = 0.1*7.5;
    driver.der_v_factor = 0.1*7.5;
  
    opt.L_kj = ...
    driver.K1 * sys.parameters.dt ...
    + driver.K2 * (rel.der_r / driver.der_r_factor)^2 ...
    + driver.K3 * (rel.der_v / driver.der_v_factor)^2;
    %
    % Continuous time dynamics
    model.opt.f    = Function('f', {sys.states.X, sys.inputs.U, sys.parameters.Z, sys.gravity.g3D}, {rel.X_dot, opt.L_kj});
    model.opt.Eq   = Function('Eq',{sys.states.X, sys.inputs.U, sys.parameters.Z, sys.gravity.g3D}, {rel.eq}, {'X','U','Z','sys.gravity.g3D'}, {'eq'});
    %Function to compute 3D gravity
    Rotsg = SX.sym('Rotsg',3,3);
    g3D_sym   = [cos(sys.states.psi), sin(sys.states.psi), 0; -sin(sys.states.psi), cos(sys.states.psi), 0; 0, 0, 1]*Rotsg*[0,0,model.global.g]';
    model.opt.fun_g3D = Function('g3D',{sys.states.psi_norm, Rotsg}, {g3D_sym}, {'sys.states.psi_norm','Rotsg'}, {'g3D_sym'});

    
    

    nx = length(sys.states.X);
    nu = length(sys.inputs.U);
    nz = length(sys.parameters.Z);

    X_scale = sys.states.X_scale;
    U_scale = [1;1;sys.inputs.delta_v_scaling];
    Z_scale = [sys.parameters.Fxf_scaling; sys.parameters.Fxr_scaling; sys.parameters.Fyf_scaling;...
               sys.parameters.Fyr_scaling;sys.parameters.Fzf_scaling;sys.parameters.Fzr_scaling;...
               sys.parameters.ax_scaling;sys.parameters.ay_scaling;1;1;sys.parameters.gamma1_scaling;...
               sys.parameters.gamma2_scaling;sys.parameters.ep_scaling;sys.parameters.dt_scaling];
  
    
    car = model;
    car.data.sys = sys;
    car.data.rel = rel;
    car.data.X_scale = X_scale;
    car.data.U_scale = U_scale;
    car.data.Z_scale = Z_scale;
    car.nx = nx;
    car.nu = nu;
    car.nz = nz;
end