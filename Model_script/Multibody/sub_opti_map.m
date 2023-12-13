function [problema, numericalData, scale]  = sub_opti_map(alfarange, pista, o, id, ID, problem_number, init_subrange, alpha_vec, opts, car)
    
% in this version we overlap only one set of states at the head and the
% tail (no controls, no algebraic variables)

    %% Discretization
    pre = diff(alfarange);
    index = find(pre<0);
    if ~isempty(index)
        pre(index) = 1-alfarange(index);
    end
    dalfa_vec = pre;
    colloc_type = car.data.colloc_type;
    d = car.data.d; 
    tau_root = [0 collocation_points(d, colloc_type)];

    %%%dalfa_vec = diff(alfarange);%round(alfarange(2) - alfarange(1),15);
    N = length(alfarange)-1;
    d = car.data.d; 
    import casadi.*    
            
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
    
    % Compute the initial state and initial guess solution
    w0 = warm_start(car,pista,d,N,index_start,index_end);

    % Set the initial state
    if problem_number == 1
        car.Xb.x0 = w0(1:car.nx);
        car.Xb.x0_lb = car.Xb.x0;
        car.Xb.x0_ub = car.Xb.x0;
    else
        car.Xb.x0 = [];
    end


    % Define numerical data for variable constraints
    Pdata.Pg = car.data.track.d_gsR_gs(:,index_start+1:index_end);
    Pdata.Pgb = car.data.track.width(:,index_start+1:index_end);
    Npg = size(Pdata.Pg,1);
    Npgb = size(Pdata.Pgb,1);
    np = 6;
    % Setup NLP
    pb = nlp_casadi_2(car.nx,car.nu,car.nz,car.nuc,car.nzc,car.nwg,N,d,car.Xb,car.Ub,car.Zb,[],[],[],Npg,Npgb,Pdata,8,np,'SXMX',colloc_type);
    
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
    % Set the initial guess solution
    pb.w0 = w0;

    % Fetch the system functions
    [~,dyn_fun,out_fun] = car.dynamic_model();

    % Dynamics constraint at collocation points
    for i = 1:d
        dx = dyn_fun(pb.xc(:,i),pb.u,pb.pg,pb.z(4:15))*pb.z(3)*car.data.Z_scale(3);
        pb.append_g(dx-pb.xp(:,i),zeros(car.nx,1),zeros(car.nx,1));
    end

    % Continuity of the solution at step endpoint
    pb.append_g(pb.xc_end-pb.x,zeros(car.nx,1),zeros(car.nx,1));

    % Value of the sparsifying variables
    v_out = out_fun(pb.x,pb.u,pb.pg);
    pb.append_g(v_out-pb.z(4:15),zeros(12,1),zeros(12,1));

    % Stay-on-track constraint
    pb.append_g(pb.x(15:17)-(pb.pg(1:3)+car.data.Z_scale(1)*pb.z(1)*pb.pg(7:9)+car.data.Z_scale(2)*pb.z(2)*pb.pg(10:12))./car.data.X_scale(15:17),zeros(3,1),zeros(3,1));
    half_track = (car.data.t1+car.data.t2)/2;
    pb.append_g(pb.z(1),-(pb.pgb(2)-half_track)/car.data.Z_scale(1),(pb.pgb(1)-half_track)/car.data.Z_scale(1));

    % Max power
    % pb.append_g(self.data.u_scale(1)*pb.u(1)*self.data.x_scale(1)*pb.x(1)/self.data.tire_radius{2}/self.data.Pa_max,0,1);
    pb.append_g(car.data.U_scale(1)*pb.u(1)*0.5*(car.data.X_scale(13)*pb.x(13)+car.data.X_scale(14)*pb.x(14))/car.data.Pa_max,0,1);

    % Complementarity constraint on input torque
    pb.append_g(pb.u(1)*pb.u(2),-5e-4,5e-4);

    % Define the cost functional
    pb.Jj = (pb.z(3)*car.data.Z_scale(3))^2 ...
        +0.001*(pb.x(2)-pb.x_1(2))^2 ... % NOTE Weight on the drift velocity variation.
        +0.001*(pb.x(6)-pb.x_1(6))^2 ... % NOTE Weight on the yaw rate variation.
        +0.010*(pb.u(1)-pb.u_1(1))^2 ... % NOTE Weight on driving torque.
        +0.010*((pb.u(2)-pb.u_1(2))*car.data.U_scale(2)/car.data.U_scale(1))^2 ... % NOTE Weight on braking torque.
        +1.000*(pb.u(3)-pb.u_1(3))^2; % NOTE Weight on steering.

    % Build map function
    pb.build_map;

    % Build g
    pb.build_g;

    % Build J
    pb.build_J;

    % Configure solver options
    opts.ipopt.linear_solver = 'ma57';
    opts.ipopt.mu_init = 1e-3;
    % opts.ipopt.alpha_for_y = 'min-dual-infeas';
    opts.ipopt.mu_strategy = 'adaptive';
    % opts.ipopt.nlp_scaling_method = 'none';
    opts.ipopt.max_iter = 1000;
    opts.ipopt.print_level = 5;
    opts.ipopt.tol = 1e-6;
    % opts.warm_start_init_point = 'yes';
    %opts.iteration_callback = MyCallback('fc',size(pb.w,1),size(pb.g,1),size(pb.p,1),pb.nx,pb.nu,pb.nz,pb.ng,N,d,self.data.x_scale,self.data.u_scale,self.data.z_scale,track);
           
                
    % NLP_casadi here is customized. After build_map it converts the
    % parameter in MX if SXMX structur is used
    if strcmp(pb.sym_type, 'SXMX')
        activation = pb.p(1);
        activation_comp = pb.p(4);
        activation_opt = pb.p(5);
        activation_start = pb.p(6);
        RHO_head = pb.p(2);
        RHO_tail = pb.p(3);
    end
    
    % Add cost function term for consensus head
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
    
    % Add cost function term for consensus tail
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
%     
%     %Solve with IPOPT
    solver = nlpsol('solver', 'ipopt', prob, opts);
    
    % Solve MLTP
%     pb.build_solver(opts,'NLP');
%     pb.solution(); 
        
    
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
    numericalData.d = d;
    
end