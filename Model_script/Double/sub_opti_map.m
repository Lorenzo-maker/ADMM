function [problema, numericalData, scale]  = sub_opti_map(alfarange, pista, o, id, ID, problem_number, init_subrange, alpha_vec, opts, car)
    
% in this version we overlap only one set of states at the head and the
% tail (no controls, no algebraic variables)

    %% Discretization
    dalfa = diff(alfarange);
    N = length(alfarange)-1;
    colloc_type = car.data.colloc_type;
    d = car.data.d;
    %% Other script
    import casadi.*
    %lap = length(find(diff(alpha_vec) < 0)) + 1;
    %car_parameters_ocp;
    %car = vehicle_casadi('double-track-full', data);
    %car = car_extra;
%     alfa_grid = alfarange;
%     
%     tau_root = [0 collocation_points(d, colloc_type)];
%     alfa_grid_colloc = zeros(d+1, N+1);
%     alfa_grid_colloc(1,:) = alfa_grid';
%     for j = 1:N
%         for i = 2:length(tau_root)
%             alfa_grid_colloc(i,j) = round(alfa_grid(j) + dalfa(j)*tau_root(i),10);
%         end
%         alfa_grid_colloc(i,N+1) = round(alfa_grid(end),10);
%     end
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
    
    pos_grid = car.data.track.pos_grid(:,index_start:index_end);%full(pista.fun_pos(alfa_grid));
    ns_grid  = car.data.track.ns_grid(:,index_start:index_end);%full(pista.fun_vh(alfa_grid));
    kp_grid  = car.data.track.kp_grid(:,index_start:index_end);%full(pista.fun_kp(alfa_grid));
    ts_grid = car.data.track.ts_grid(:,index_start:index_end);%full(pista.fun_ts(alfa_grid));
    psi_grid = car.data.track.psi_grid(:,index_start:index_end);%full(pista.fun_ts(alfa_grid));
%     ts_grid(3,:) = 0;
%     ts_grid = ts_grid./vecnorm(ts_grid);
%     if length(ID.H) > 1
%         if isempty(ID.H{problem_number}) && isempty(ID.h{problem_number})
%             out = [ID.T{problem_number},ID.t{problem_number}];
%             psi_grid = psi_grid(1:max(out));
%         elseif isempty(ID.T{problem_number}) && isempty(ID.t{problem_number})
%             in = [ID.H{problem_number},ID.h{problem_number}];
%             psi_grid = psi_grid(min(in):end);
%         else
%             in = [ID.H{problem_number},ID.h{problem_number}];
%             out = [ID.T{problem_number},ID.t{problem_number}];
%             psi_grid = psi_grid(min(in):max(out));
%         end
%     end
    RVSG = car.data.track.RVSG(index_start:index_end);
    RVSG_colloc = car.data.track.RVSG_colloc(:,index_start:index_end-1);
%     RVSG = cell(1,N+1);
%     RVSG_colloc = cell(d,N+1);
%     for i = 1:N+1
%         RVSG{i} = [cos(-psi_grid(i)), sin(-psi_grid(i)), 0; -sin(-psi_grid(i)), cos(-psi_grid(i)), 0; 0, 0, 1]*(full(pista.fun_Rgs(alfa_grid(i))))';
%         for j = 1:d
%             RVSG_colloc{j, i} = [cos(-psi_grid(i)), sin(-psi_grid(i)), 0; -sin(-psi_grid(i)), cos(-psi_grid(i)), 0; 0, 0, 1]*(full(pista.fun_Rgs(alfa_grid_colloc(j+1,i))))';
%         end
%     end
%     RVSG_colloc = RVSG_colloc(:, 1:end-1); % N dimension
    %psi_grid = atan_track(ts_grid, 'clockwise');
    %% Construct NLP
    if problem_number > 1
        %%%%% In opti singolo ef = 0 %%%%%%%
        car.Xb.x0 = [];        
    else
        car.Xb.x0 = [car.data.Vi;0;0;pos_grid(1,1);pos_grid(2,1);psi_grid(1);repmat(car.data.Vi/car.data.rw,2,1); repmat(car.data.Vi*(1 - 0)/car.data.rw,2,1)]./car.data.X_scale;
        car.Xb.x0_lb = [car.data.Vi;0;0;pos_grid(1,1);pos_grid(2,1);psi_grid(1);repmat(car.data.Vi/car.data.rw,2,1); repmat(car.data.Vi*(1 - 0.05)/car.data.rw,2,1)]./car.data.X_scale;
        car.Xb.x0_ub = [car.data.Vi;0;0;pos_grid(1,1);pos_grid(2,1);psi_grid(1);repmat(car.data.Vi/car.data.rw,2,1); repmat(car.data.Vi*(1 + 0.05)/car.data.rw,2,1)]./car.data.X_scale;
    end
    x00 = [car.data.Vi;0;0;pos_grid(1,1);pos_grid(2,1);psi_grid(1);repmat(car.data.Vi/car.data.rw,2,1); repmat(car.data.Vi*(1 - 0)/car.data.rw,2,1)]./car.data.X_scale;
    np = 6;
    % Numerical data for variable constraints
    Pdata.Pgb = [-pos_grid(4,2:end); pos_grid(4,2:end)];
    Rvsg = RVSG;
    Rvsg = horzcat(Rvsg{:});
    Rvsg = Rvsg(:);
    Rvsg = reshape(Rvsg, 9, []);
    Rvsg_colloc = horzcat(RVSG_colloc{:});
    Rvsg_colloc = Rvsg_colloc(:);
    Rvsg_colloc = reshape(Rvsg_colloc, 9*d, []);
    Pdata.Pg = [pos_grid(1:2,2:end); ns_grid(1:2, 2:end); Rvsg(:,2:end); Rvsg_colloc];

    % NLP
    pb = nlp_casadi_2(car.nx, car.nu, car.nz, car.nuc, car.nzc, car.nwg, N, d, car.Xb, car.Ub, car.Zb, [], [], [], 4 + 9 + 9*d, 2, Pdata, 8, np, 'SXMX', colloc_type);
    pb.w0 = [];
    for k = 1:N
        initial_guess_k;    
        pb.w0 = [pb.w0;w0_k];
    end
    pb.w0 = [x00; pb.w0];
    
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
    pb.append_g(pb.x_1(9) - pb.x_1(10), 0, 0, 'start');
    RvsgKj = reshape(pb.pg(5:end),9, d+1);
    RvsgK = reshape(RvsgKj(:, 1),3,3);
    for i = 1:d   
       Rvsgj = reshape(RvsgKj(:, i+1),3,3);
       g3DKj = car.fun_g3D(pb.xc(6,i), Rvsgj); 
       pb.append_g(car.F(pb.xc(:,i), pb.u, pb.z, g3DKj)*pb.z(end)*car.data.Z_scale(end) - pb.xp(:,i).*car.data.X_scale, zeros(car.nx,1), zeros(car.nx,1)) ;
    end
    g3DK = car.fun_g3D(pb.x(6), RvsgK); 
    xdot = car.F(pb.x, pb.u, pb.z, g3DK);
    pb.Jj = (pb.z(end)*car.data.Z_scale(end))^2 + 10*((pb.u(3) - pb.u_1(3))*car.data.tau)^2 + 0.0001*xdot(2)^2 + 0.01*xdot(3)^2 + 0*0.0000001*(xdot(end-4)^2 + xdot(end-3)^2 + xdot(end-2)^2 + xdot(end-1)^2);
    % Continuity equation for generic interval
    pb.append_g(pb.xc_end - pb.x , zeros(car.nx,1), zeros(car.nx,1));
    % Stay on track constraint
    pb.append_g(pb.x(4:5) - (pb.pg(1:2) + pb.pg(3:4)*pb.z(end-1)*car.data.Z_scale(end-1))./car.data.X_scale(4:5), zeros(2,1), zeros(2,1));
    pb.append_g(pb.z(end-1), pb.pgb(1)/car.data.Z_scale(end-1), pb.pgb(2)/car.data.Z_scale(end-1));    
    % Algebraic equations constraints
    pb.append_g(car.Eq(pb.x, pb.u, pb.z, g3DK), zeros(car.neq,1), zeros(car.neq,1));
    % Steering constraint
    % pb.append_g(pb.u(3) - pb.u_1(3), -pi/data.U_scale(3), pi/data.U_scale(3));
    % Complementary constraint
    pb.append_g(pb.u(1)*pb.u(2), -1e-6, 1e-6);
    % PowerTrain
    pb.append_g(((pb.u(1)*car.data.U_scale(1)*pb.x(1)*car.data.X_scale(1)/car.data.rw) + (pb.u(2)*car.data.U_scale(2)*pb.x(1)*car.data.X_scale(1)/car.data.rw))/car.data.P_scale, car.data.Pmin/car.data.P_scale, car.data.Pmax/car.data.P_scale);
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
    
    %Solve with IPOPT
    solver = nlpsol('solver', 'ipopt', prob, opts);
    



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