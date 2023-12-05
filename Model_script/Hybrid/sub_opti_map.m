function [problema, numericalData, scale]  = sub_opti_map(alfarange, pista, o, id, ID, problem_number, init_subrange, alpha_vec, opts, car)
    
% in this version we overlap only one set of states at the head and the
% tail (no controls, no algebraic variables)

    %% Discretization
    N = length(alfarange)-1;
    
    %% Other script
    import casadi.*
    lap = length(find(diff(alpha_vec) < 0)) + 1;
    Nsteps = length(alpha_vec)-1; % total step of the whole grid
    d = car.opt.d_colloc;
    colloc_type = car.opt.colloc_type;
    %% Model Parameters (Build Model Structure)
    %car_parameters_ocp
    %[sys, rel] = fnc_vehicleModel(car, pista);
    sys = car.data.sys;
    rel = car.data.rel;

    
%     driver.K1 = 0.95;
%     driver.K2 = 0.025;
%     driver.K3 = 0.025;
%     driver.der_r_factor = 0.1*7.5;
%     driver.der_v_factor = 0.1*7.5;
    
    
    alfa_grid = alfarange;
    if length(alfarange) == ID.end
        index_start = ID.start;
        index_end = ID.end;
    else
        if problem_number > 1 && problem_number < length(ID.H)
            index_start = ID.H{problem_number}(1);
            index_end = ID.T{problem_number}(end);
        elseif problem_number == 1
            index_start = ID.start;
            index_end = ID.T{problem_number}(end);
        elseif problem_number == length(ID.H)
            index_start = ID.H{problem_number}(1);
            index_end = ID.end;
        end
    end
    pos_grid = car.track.pos_grid(:,index_start:index_end);%full(pista.fun_pos(alfa_grid));
    ns_grid  = car.track.ns_grid(:,index_start:index_end);%full(pista.fun_vh(alfa_grid));
    kp_grid  = car.track.kp_grid(:,index_start:index_end);%full(pista.fun_kp(alfa_grid));
    ts_grid = car.track.ts_grid(:,index_start:index_end);%full(pista.fun_ts(alfa_grid));
    psi_grid = car.track.psi_grid(:,index_start:index_end);%full(pista.fun_ts(alfa_grid));    
    RVSG = car.track.RVSG(index_start:index_end);
    
%     % Objective term
%     opt.L_kj = ...
%         driver.K1 * sys.parameters.dt ...
%         + driver.K2 * (rel.der_r / driver.der_r_factor)^2 ...
%         + driver.K3 * (rel.der_v / driver.der_v_factor)^2;
%     %
%     % Continuous time dynamics
%     opt.f    = Function('f', {sys.states.X, sys.inputs.U, sys.parameters.Z, sys.gravity.g3D}, {rel.X_dot, opt.L_kj});
%     opt.Eq   = Function('Eq',{sys.states.X, sys.inputs.U, sys.parameters.Z, sys.gravity.g3D}, {rel.eq}, {'X','U','Z','sys.gravity.g3D'}, {'eq'});
%     %Function to compute 3D gravity
%     Rotsg = SX.sym('Rotsg',3,3);
%     g3D_sym   = [cos(sys.states.psi), sin(sys.states.psi), 0; -sin(sys.states.psi), cos(sys.states.psi), 0; 0, 0, 1]*Rotsg*[0,0,model.global.g]';
%     fun_g3D = Function('g3D',{sys.states.psi_norm, Rotsg}, {g3D_sym}, {'sys.states.psi_norm','Rotsg'}, {'g3D_sym'});
    
    opt.f = car.opt.f;
    opt.Eq = car.opt.Eq;
    fun_g3D = car.opt.fun_g3D;
    model = car;
    %% Construct NLP
    if problem_number > 1
        %%%%% In opti singolo ef = 0 %%%%%%%
        opt.Xb.x0 = [];        
    else
        opt.Xb.x0 = sys.states.X_0;
    end
    opt.Xb.lb = sys.states.X_lb;
    opt.Xb.ub = sys.states.X_ub;
    opt.Ub.lb = sys.inputs.U_lb;
    opt.Ub.ub = sys.inputs.U_ub;
    opt.Zb.lb = sys.parameters.Z_lb;
    opt.Zb.ub = sys.parameters.Z_ub;
    opt.Npg = 3*3 + 2 + 2;
    opt.Npgb = 2;
    opt.Nthread = 8;
    opt.np = 6;
    % Numerical data for variable constraints
%     RVSG = cell(1,N+1);
%     for i = 1:N+1
%         RVSG{i} = [cos(-psi_grid(i)), sin(-psi_grid(i)), 0; -sin(-psi_grid(i)), cos(-psi_grid(i)), 0; 0, 0, 1]*(full(pista.fun_Rgs(alfa_grid(i))))';
%     end    
    for k = 1:N
        opt.data_g.Pg(:,k) = [RVSG{k+1}(:); pos_grid(1:2,k+1); ns_grid(1:2,k+1)];
    end
    opt.data_g.Pgb = [-(pos_grid(4,2:end)-model.opt.ep_margin);pos_grid(4,2:end)-model.opt.ep_margin];    

    % NLP
    pb = nlp_casadi_2(sys.states.nx, sys.inputs.nu, sys.parameters.np, 0, 0, 0, N, d, opt.Xb, opt.Ub, opt.Zb, [], [], [], opt.Npg, opt.Npgb, opt.data_g, opt.Nthread, opt.np, 'SXMX', colloc_type);
    for k = 1:N
        initial_guess
        if k == 1
            if isempty(opt.Xb.x0)
                pb.w0 = [pb.w0; sys.states.X_0; sys.inputs.U_0; sys.parameters.Z_0; repmat(sys.states.X_0, model.opt.d_colloc+1, 1)];
            else
                pb.w0 = [pb.w0; opt.Xb.x0; sys.inputs.U_0; sys.parameters.Z_0; repmat(sys.states.X_0, model.opt.d_colloc+1, 1)];
            end
        else
            pb.w0 = [pb.w0; sys.inputs.U_0; sys.parameters.Z_0; repmat(sys.states.X_0, model.opt.d_colloc+1, 1)]; 
        end
    end
    
    if strcmp(pb.sym_type, 'SXMX')
        sym_type = 'MX';
    else
        sym_type = pb.sym_type;
    end
    activation = pb.p(1);
    activation_comp = pb.p(4);
    activation_opt = pb.p(5);
    activation_start = pb.p(6);
    RHO_head = pb.p(2);
    RHO_tail = pb.p(3); 
    
    % Collocation
    % NLP
    Rgs = reshape(pb.pg(1:9),3,3);
    for i = 1:model.opt.d_colloc
        g3DKj = fun_g3D(pb.xc(6,i),Rgs);
        [opt.fr, opt.qr] = opt.f(pb.xc(:,i), pb.u, pb.z, g3DKj);
        pb.append_g(opt.fr*pb.z(end)*sys.parameters.dt_scaling - pb.xp(:,i).*sys.states.X_scale, zeros(pb.nx,1), zeros(pb.nx,1)) ;
        pb.Jj = [pb.Jj; opt.qr * pb.z(end)*sys.parameters.dt_scaling];
    end
    pb.Jj = pb.Jj*activation_opt + activation_start*(pb.z(11)*sys.parameters.ep_scaling)^2;
    % Evaluate cost function for generic interval
    pb.cost;
    % Driver
    pb.append_g((pb.u(3) - pb.u_1(3))*sys.inputs.delta_v_scaling/(pb.z(end)*sys.parameters.dt_scaling + 1e-5), -pi, pi);
    % Continuity equation for generic interval
    pb.append_g(pb.xc_end - pb.x, zeros(pb.nx,1), zeros(pb.nx,1));
    % Throttle/Braking
    pb.append_g(pb.u(1)*pb.u(2), 0, 1e-5);
    % beta_r
    pb.append_g(pb.z(9) - pb.u(2), -10, 0);
    % gamma_r
    pb.append_g(pb.z(10) - pb.u(1), -10, 0);
    %adherence_f
    pb.append_g((pb.z(1) * sys.parameters.Fxf_scaling / ((model.vehicle.tire.mux_0(1) + model.vehicle.tire.k_mux(1) * pb.z(5) * sys.parameters.Fzf_scaling) * pb.z(5) * sys.parameters.Fzf_scaling) ...
            * 0.5 * (model.opt.gripxa_max(1)^-1 + model.opt.gripxb_max(1)^-1 + if_else_smooth(pb.z(1),0,1,-1) * (model.opt.gripxa_max(1)^-1 - model.opt.gripxb_max(1)^-1)))^2 ...
            + (pb.z(3) * sys.parameters.Fyf_scaling / (model.opt.gripy_max(1) * (model.vehicle.tire.muy_0(1) + model.vehicle.tire.k_muy(1) * pb.z(5) * sys.parameters.Fzf_scaling) * pb.z(5) * sys.parameters.Fzf_scaling))^2, 0, 1);
    %adherence_r
    pb.append_g((pb.z(2) * sys.parameters.Fxr_scaling / ((model.vehicle.tire.mux_0(2) + model.vehicle.tire.k_mux(2) * pb.z(6) * sys.parameters.Fzr_scaling) * pb.z(6) * sys.parameters.Fzr_scaling) ...
            * 0.5 * (model.opt.gripxa_max(2)^-1 + model.opt.gripxb_max(2)^-1 + if_else_smooth(pb.z(2),0,1,-1) * (model.opt.gripxa_max(2)^-1 - model.opt.gripxb_max(2)^-1)))^2 ...
            + (pb.z(4) * sys.parameters.Fyr_scaling / (model.opt.gripy_max(2) * (model.vehicle.tire.muy_0(2) + model.vehicle.tire.k_muy(2) * pb.z(6) * sys.parameters.Fzr_scaling) * pb.z(6) * sys.parameters.Fzr_scaling))^2, 0, 1);
    %ep
    pb.append_g(pb.z(11)*sys.parameters.ep_scaling, pb.pgb(1), pb.pgb(2));
    %track
    pos = pb.pg(10:11);
    n = pb.pg(12:13);
    pb.append_g(pb.x(4:5) - (pos + n*pb.z(11)*sys.parameters.ep_scaling)./sys.states.X_scale(4:5), zeros(2,1), zeros(2,1))
    %Eq
    g3DK = fun_g3D(pb.x(6),Rgs);
    pb.append_g(vertcat(opt.Eq(pb.x, pb.u, pb.z, g3DK)), zeros(rel.neq,1), zeros(rel.neq,1));
    %terminal
    if problem_number == length(ID.H)
        pb.append_g((model.vehicle.bat.SOC_min + pb.x(7) * sys.states.SOC_scaling) - model.opt.SOC_final, 0, 1, 'end');
    end
    % % Cyclic
    % pb.append_g([pb.x0(1:3);pb.z0(11)]-[pb.x_end(1:3);pb.z_end(11)], [-1;zeros(3,1)], [1;zeros(3,1)], 'cyclic');
    % Build map function 
    pb.build_map;
    % Build g 
    pb.build_g;
    % Build J
    pb.build_J;

    
    if strcmp(pb.sym_type, 'SXMX')
        activation = pb.p(1);
        activation_comp = pb.p(4);
        activation_opt = pb.p(5);
        activation_start = pb.p(6);
        RHO_head = pb.p(2);
        RHO_tail = pb.p(3);
    end
    
    
    Jc = 0;
    Zcons = [];
    Ycons = [];
    if ~isempty(id.h{problem_number})
        if o == 0
            k_head = id.h{problem_number};
        else
            k_head = id.h{problem_number}(1:end-1);
        end
        % consensus term        
        Zcons_h = casadi.(sym_type).sym(['Zcons' num2str(k_head)], o*(car.nx + car.nu + car.nz) + car.nx);
        Zcons = [Zcons; Zcons_h];        
        Ycons_h = casadi.(sym_type).sym(['Ycons' num2str(k_head)], o*(car.nx + car.nu + car.nz) + car.nx);
        Ycons = [Ycons; Ycons_h];
        if o == 0
            x_h = pb.X_1(:,k_head);
        else
            x_h = [pb.X_1(:,k_head); pb.U(:,k_head); pb.Z(:,k_head)];
            x_h = [x_h(:); pb.X(:,id.h{problem_number}(end-1))];
        end
        cdiff = (x_h - Zcons_h);
        Jc = Jc + Ycons_h'*cdiff + 0.5*RHO_head*(cdiff.')*(cdiff);
    end
    
    if ~isempty(id.t{problem_number})
        if o == 0
            k_tail = id.t{problem_number};
        else            
            k_tail = id.t{problem_number}(1:end-1);
        end
        % consensus term        
        Zcons_t = casadi.(sym_type).sym(['Zcons' num2str(k_tail)], o*(car.nx + car.nu + car.nz) + car.nx);
        Zcons = [Zcons; Zcons_t];
        Ycons_t = casadi.(sym_type).sym(['Ycons' num2str(k_tail)], o*(car.nx + car.nu + car.nz) + car.nx);
        Ycons = [Ycons; Ycons_t];        
        if o == 0
            x_t = pb.X(:,k_tail-1);
        else            
            x_t = [pb.X_1(:,k_tail); pb.U(:,k_tail); pb.Z(:,k_tail)];
            x_t = [x_t(:); pb.X(:,id.t{problem_number}(end-1))];
        end
        cdiff = (x_t - Zcons_t);
        Jc = Jc + Ycons_t'*cdiff + 0.5*RHO_tail*(cdiff.')*(cdiff);
    end
    pb.J = pb.J + activation.*(Jc);
    pb.p = [Zcons; Ycons; pb.p];
    
    
    % Create an NLP solver 
    prob = struct('f', pb.J, 'x', pb.w, 'g', pb.g, 'p', pb.p);
    
    options = struct( ...
    'ipopt', ...
    struct( ...
    'max_iter', 5000, ... 
    'linear_solver', 'ma57', ...
    'print_level', 5, ... % 
    'mu_strategy', 'adaptive', ...
    'ma57_pre_alloc', 2.5, ...
    'nlp_scaling_method', 'none', ...
    'tol', 1e-6), ...
    'print_time', true);

    %Solve with IPOPT
    solver = nlpsol('solver', 'ipopt', prob, options);
    



    %Solve with WORHP
%     opts = struct;
%     opts.worhp.Algorithm = 1;%S
%     opts.worhp.MA97blas3 = true;%s
%     opts.worhp.MA97mf = true;%S
%     opts.worhp.MaxIter = 200;%s
% %     opts.worhp.qp.ipLsMethod = 4;%s
%     opts.worhp.MA97scaling = 0;%s
%     opts.print_time = 0;%s
%     opts.worhp.qp.ipComTol = 1e-5;
%     opts.worhp.TolFeas = 1e-5;
%     opts.worhp.TolOpti = 1e-5;
%     opts.worhp.TolComp = 1e-5;

%     
%     
%     opts.worhp.NLPprint = -1;
%     opts.worhp.qp.lsScale = false;%s
%     % opts.worhp.qp.method = 12;
%     opts.worhp.qp.ipBarrier = 2;%s
%     opts.worhp.ScaledQP = false;%s
%     opts.worhp.ScaledObj = false;%s
    % opts.calc_lam_p = false;
    % opts.worhp.IP_BarrierInit = 1e-3;
    % opts.worhp.Crossover = 2;
    % opts.worhp.CrossoverTol  = 0.1;
    % opts.worhp.FeasibleDual = true;

%     solver = nlpsol('solver', 'worhp', prob, opts);
    
%     batFile = sprintf('%s_%i.bat', 'compilation', problem_number);
%     c_filename = sprintf('%s_%i.c', 'sub_problem', problem_number);
%     solver.generate_dependencies(c_filename);
%     diskName = pwd;
%     diskName = diskName(1:2);
%     fid = fopen(batFile, 'w');
%     fprintf(fid, '%s\n', diskName);  
%     fprintf(fid, 'cd %s\\%s\n', pwd, 'Compilati');
%     fprintf(fid, 'cl /LD %s\n', c_filename);
%     fclose(fid);
%     
%     movefile(c_filename, 'Compilati');
%     movefile(batFile, 'Compilati');
    
    problema = solver;
    numericalData.w0 = pb.w0;
    numericalData.initialCons = zeros(length(Zcons), 1);
    numericalData.lbw = pb.lbw;
    numericalData.ubw = pb.ubw;
    numericalData.lbg = pb.lbg;
    numericalData.ubg = pb.ubg;
    numericalData.nx = car.nx;
    numericalData.nu = car.nu;
    numericalData.nz = car.nz;
    numericalData.opts = opts;
%     numericalData.dll_filename = sprintf('%s_%i.dll', 'Compilati\sub_problem', problem_number);
%     numericalData.bat_filename = sprintf('%s_%i.bat', 'Compilati\compilation', problem_number);  
    numericalData.d = d;
    
end