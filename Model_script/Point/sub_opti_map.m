function [problema, numericalData, scale]  = sub_opti_map(alfarange, pista, o, id, ID, problem_number, init_subrange, alpha_vec, opts, car)
    
% in this version we overlap only one set of states at the head and the
% tail (no controls, no algebraic variables)

    %% Discretization
    dalfa = round(alfarange(2) - alfarange(1),15);
    N = length(alfarange)-1;
    colloc_type = car.data.colloc_type;
    d = car.data.d;
    %% Other script
    import casadi.*
    lap = length(find(diff(alpha_vec) < 0)) + 1;
%     car_parameters_ocp;
%     car = vehicle_casadi('point-mass', data);
%     alfa_grid = alfarange;
%     pos_grid = full(pista.fun_pos(alfa_grid));
%     ns_grid  = full(pista.fun_vh(alfa_grid));
%     kp_grid  = full(pista.fun_kp(alfa_grid));
%     ts_grid = full(pista.fun_ts(alfa_grid));
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
    %psi_grid = atan_track(ts_grid, 'clockwise');
    
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
    %% Construct NLP
    if problem_number > 1
        %%%%% In opti singolo ef = 0 %%%%%%%
        car.Xb.x0 = [];        
    else
        car.Xb.x0 = [car.data.Vi;0;0;pos_grid(1,1);pos_grid(2,1);psi_grid(1)]./car.data.X_scale;
        car.Xb.x0_lb = [car.data.Vi;0;0;pos_grid(1,1);pos_grid(2,1);psi_grid(1)]./car.data.X_scale;
        car.Xb.x0_ub = [car.data.Vi;0;0;pos_grid(1,1);pos_grid(2,1);psi_grid(1)]./car.data.X_scale;
    end
    np = 6;
    % Numerical data for variable constraints
    Pdata.Pgb = [-pos_grid(4,2:end); pos_grid(4,2:end)];
    Pdata.Pg = [pos_grid(1:2,2:end); ns_grid(1:2, 2:end)];

    % NLP
    pb = nlp_casadi_2(car.nx, car.nu, car.nz, car.nuc, car.nzc, car.nwg, N, d, car.Xb, car.Ub, car.Zb, [], [], [], 4, 2, Pdata, 8, np, 'SXMX', colloc_type);
    pb.w0 = [];
    for k = 1:N
        initial_guess_k;    
        pb.w0 = [pb.w0;w0_k];
    end
    pb.w0 = [[car.data.Vi;0;0;pos_grid(1,1);pos_grid(2,1);psi_grid(1)]./car.data.X_scale; pb.w0];
    
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
    for i = 1:d
       pb.append_g(car.F(pb.xc(:,i), pb.u, pb.z)*pb.z(5)*car.data.Z_scale(5) - pb.xp(:,i).*car.data.X_scale, zeros(car.nx,1), zeros(car.nx,1)) ;
    end
    xdot = car.F(pb.x, pb.u, pb.z);
    pb.Jj = activation_opt*((pb.z(5)*car.data.Z_scale(5))^2 + 0.1*(pb.u(2)*car.data.U_scale(2))^2 + 0.01*xdot(2)^2 + 1*(pb.z_1(2)-pb.z(2))^2);
    pb.Jj = pb.Jj + activation_start*((pb.x(1)*car.data.X_scale(1) - 20)^2 + (pb.z(4)*car.data.Z_scale(4))^2);
    pb.append_g(pb.xc_end - pb.x, zeros(car.nx,1), zeros(car.nx, 1));

    % pb.append_g((pb.x - pb.x_1).*data.X_scale - xdot.*pb.z(5)*data.Z_scale(5), zeros(car.nx,1), zeros(car.nx,1))
    % Stay on track constraint
    pb.append_g(pb.x(4:5) - (pb.pg(1:2) + pb.pg(3:4)*pb.z(4)*car.data.Z_scale(4))./car.data.X_scale(4:5), zeros(2,1), zeros(2,1));
    pb.append_g(pb.z(4), pb.pgb(1)/car.data.Z_scale(4), pb.pgb(2)/car.data.Z_scale(4));    
    % Adherence constraint
    pb.append_g((pb.z(2)*car.data.Z_scale(2)/(car.data.mu_y*pb.z(1)*car.data.Z_scale(1)))^2 + (pb.u(1)*car.data.U_scale(1)/(car.data.mu_x*pb.z(1)*car.data.Z_scale(1)))^2,0,1); 
    % Algebraic equations constraints
    pb.append_g(car.Eq(pb.x, pb.u, pb.z), zeros(car.neq,1), zeros(car.neq,1));
   
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