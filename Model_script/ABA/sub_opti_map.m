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
    tau_colloc = tau_root(2:end);

    %%%dalfa_vec = diff(alfarange);%round(alfarange(2) - alfarange(1),15);
    N = length(alfarange)-1;
    colloc_type = car.data.colloc_type;
    d = car.data.d; 
    import casadi.*
    %% Unpack car data 
    nx = car.data.nx; nu = car.data.nu; nz = car.data.nz;
    X_lb = car.data.X_lb; X_ub = car.data.X_ub;
    U_lb = car.data.U_lb; U_ub = car.data.U_ub;
    Z_lb = car.data.Z_lb; Z_ub = car.data.Z_ub;
    X_scale = car.data.X_scale;
    U_scale = car.data.U_scale;
    Z_scale = car.data.Z_scale;
    ep_scale = car.data.ep_scale;
    Fy_scale = car.data.Fy_scale;
    Mz_scale = car.data.Mz_scale;
    Fx_scale = car.data.Fx_scale;
    Vi = car.data.Vi;
    Vmax = car.data.Vmax;
    kb = car.data.kb;
    t = car.data.t;
    rho = car.rho;
    S = car.Sf;
    cx = car.cx;
    cz1 = car.cz1;
    cz2 = car.cz2;    
    Pmin = car.data.Pmin;
    Pmax = car.data.Pmax;
    tol = car.data.tol;
    F = car.fun.F;
    Fz11tyre = car.fun.Fz11tyre;
    Fz12tyre = car.fun.Fz12tyre;
    Fz21tyre = car.fun.Fz21tyre;
    Fz22tyre = car.fun.Fz22tyre;
    Fytyre = car.fun.Fytyre;
    Fymaxij = car.fun.Fymaxij;
    Fxmaxij = car.fun.Fxmaxij;
    alfa1_fun = car.fun.alfa1_fun;
    alfa2_fun = car.fun.alfa2_fun;
    DZ_1 = car.fun.DZ_1;
    DZ_2 = car.fun.DZ_2;
    
    %% Build index of start and end subproblem (w.r.t. the whole grid)
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
    
    %% Exctract numerical track data (already stored in the model)
    g_gs_num = car.data.track.g_gs_num(:,:, index_start:index_end);
    Jac_num = car.data.track.Jac_num(:, index_start:index_end);
    Jac_der_num = car.data.track.Jac_der_num(:, index_start:index_end);
    
    dalfa_vec = car.data.track.dalfa_vec;
    g_gs_colloc = cell(length(alfarange)-1, d);
    Jac_colloc = cell(length(alfarange)-1, d);
    Jac_der_colloc = cell(length(alfarange)-1, d);
    %index_colloc = linspace(1,length(alfarange)-1,length(alfarange)-1) + (index_start-1);
    index_colloc = linspace(index_start,index_end-1, index_end-index_start);
    for k = 1:length(alfarange)-1
        for j = 1:d    
            g_gs_colloc{k,j} = car.data.track.g_gs_colloc{index_colloc(k),j};
            Jac_colloc{k,j} = car.data.track.Jac_colloc{index_colloc(k),j};
            Jac_der_colloc{k,j} = car.data.track.Jac_der_colloc{index_colloc(k),j};
        end
    end

%     g_gs_num = full(pista.fun_g_gs(alfarange));
%     g_gs_num = reshape(g_gs_num, 4, 4, []);
%     Jac_num = full(car.T_jac(alfarange));
%     Jac_der_num = chop(full(car.T_jac_der(alfarange)),1);
%     
%     %Evaluate numerical quantities
%     g_gs_colloc = cell(length(alfarange)-1, d);
%     Jac_colloc = cell(length(alfarange)-1, d);
%     Jac_der_colloc = cell(length(alfarange)-1, d);
%     alfa_colloc = zeros(length(alfarange)-1, d);
%     for k = 1:length(alfarange)-1
%         for j = 1:d
%             alfa_colloc(k,j) = min(alfarange(k) + dalfa_vec(k)*tau_colloc(j),1);
%             g_gs_colloc{k, j} = full(pista.fun_g_gs(min(alfarange(k) + dalfa_vec(k)*tau_colloc(j),1)));
%             Jac_colloc{k, j} = full(car.T_jac(min(alfarange(k) + dalfa_vec(k)*tau_colloc(j),1)));
%             Jac_der_colloc{k, j} = chop(full(car.T_jac_der(min(alfarange(k) + dalfa_vec(k)*tau_colloc(j),1))),1);
%         end
%     end


    %% Import initial guess
    if ~isempty(init_subrange)
        guess_0 = init_subrange;
        guess_0.q = guess_0.q(:, 2:end); % remove alpha from guess
        X_0 = guess_0.q(1,:)'./X_scale; 
        guess_0.q = guess_0.q;
    else
        X_0 = zeros(nx,1);
        X_0(3) = -0.022;
        X_0(6) = 6e-4;
        X_0 = X_0./X_scale;
    end
    %% Build inpuit for NLP casadi    
    % Define lower & upper bound
    Xb.lb = X_lb;
    Xb.ub = X_ub;
    if problem_number > 1
        %%%%% In opti singolo ef = 0 %%%%%%%
        Xb.x0 = [];        
    else
        Xb.x0 = X_0;
        Xb.x0_lb = ([X_lb(1); 0; X_lb(3:end)]);
        Xb.x0_ub = ([X_ub(1); 0; X_ub(3:end)]);
    end
    Ub.lb = U_lb;
    Ub.ub = U_ub;
    Zb.lb = Z_lb;
    Zb.ub = Z_ub;
    
    % Define number of symbolic and numerical parameters
    np = 6;
    Npg = (4*4 + 6 + 6 + 1)*(d+2);
    Npgb = 2;
    % Evaluating numerical data from track (maybe can be stored in car.data)
    Pdata = struct;
    Pdata.Pg = zeros((16 + 6 + 6 + 1)*(d+2), N);
    for k = 1:N
        Pdata.Pg(1:(16+6+6+1),k) = [reshape(g_gs_num(:, :, k),16,1); Jac_num(:, k); Jac_der_num(:, k); dalfa_vec(k)];
        for j = 1:d
            Pdata.Pg((16+6+6+1) + (j-1)*(16+6+6+1) + 1:(16+6+6+1) + j*(16+6+6+1),k) = [g_gs_colloc{k, j}(:); Jac_colloc{k, j}; Jac_der_colloc{k, j}; dalfa_vec(k)];
        end
        Pdata.Pg(end-(16+6+6+1)+1:end,k) = [reshape(g_gs_num(:, :, k+1),16,1); Jac_num(:, k+1); Jac_der_num(:, k+1); dalfa_vec(k)];
    end
    ep_pos = car.data.track.pos(5,index_start:index_end);%pista.fun_pos(alfarange);
    Pdata.Pgb = [(ep_pos(2:end) - t/2); ep_pos(2:end) - t/2]./ep_scale;
    
    %% Build pb object with NLP casdi
    pb = nlp_casadi_2(nx, nu, nz, 0, 0, 0, N, d, Xb, Ub, Zb, [], [], [], Npg, Npgb, Pdata, 8, np, 'SXMX', colloc_type);
    
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
    
    % Initial guess
    for k = 1:N        
        V0 = Vi;
        X_0 = zeros(11,1);
        if ~isempty(init_subrange)
            X_0(3) = mean(guess_0.q(:,3))/X_scale(3);
            X_0(6) = mean(guess_0.q(:,6))/X_scale(6);
        else
            X_0(3) = -0.022/X_scale(3);
            X_0(6) = 6e-04/X_scale(6);
        end
        Fz0 = car.M(1,1)*9.81 + 0.5*rho*(cz1+cz2)*S*V0^2;
        U_0 = [0.5*rho*cx*S*V0^2; 0; car.M(1,1)*9.81 + 0.5*rho*(cz1+cz2)*S*V0^2; 0; 0; 0; 0]./U_scale;
        Z_0 = [Fz0/4*ones(4,1); 0; 0; 0.5*rho*cx*S*V0^2; 0]./Z_scale;
        pb.w0 = [pb.w0; U_0; Z_0; repmat(X_0, d+1, 1)];
    end
    if ~isempty(init_subrange)
        pb.w0 = [guess_0.q(1,:)'./X_scale; pb.w0];
    else
        pb.w0 = [X_0./X_scale; pb.w0];
    end
    
    % Start constraint (alfa = 0)
    Vk = car.axle_twist(g_gs_num(:, :, 1), Jac_num(:, 1), Jac_der_num(:, 1), pb.x_1.*X_scale); 
    if problem_number > 1
        pb.append_g(Vk(1)/Vmax, 0, 1, 'start');
    else
        pb.append_g(Vk(1)/Vmax, Vi/Vmax, Vi/Vmax, 'start');               
    end
    pb.append_g(pb.x_1(1), -(ep_pos(1)-t/2)/ep_scale, (ep_pos(1)-t/2)/ep_scale, 'start');
     
    % Collocation
    data_num = reshape(pb.pg, 16+6+6+1, d+2);
    dalfa = data_num(end,1);
    data_num = data_num(1:end-1,:);
    for i = 1:d
        g_gs_colloc_j = reshape(data_num(1:16, i+1),4,4);
        Jac_colloc_j = (data_num(16 + 1: 16 + 6, i+1));
        Jac_der_colloc_j = (data_num(16 + 7: end, i+1));
        fj = F(g_gs_colloc_j, Jac_colloc_j, Jac_der_colloc_j, pb.xc(:,i).*X_scale, pb.u.*U_scale);
        pb.append_g((dalfa*fj(2:end) - pb.xp(:,i).*X_scale)./X_scale, zeros(nx,1), zeros(nx, 1));
    end
    
    % Continuity
    pb.append_g(pb.xc_end - pb.x, zeros(nx,1), zeros(nx, 1));
    
    % Obtaining numerical quantities for the generic step
    g_gs_num_k = reshape(data_num(1:16, d+2),4,4);
    Jac_num_k = (data_num(16 + 1: 16 + 6, d+2));
    Jac_der_num_k = (data_num(16 + 7: end, d+2));
    
    Vk = car.axle_twist(g_gs_num_k, Jac_num_k, Jac_der_num_k, pb.x.*X_scale);  
    g_gs_num_1 = reshape(data_num(1:16, 1),4,4);
    Jac_num_1 = (data_num(16 + 1: 16 + 6, 1));
    Jac_der_num_1 = (data_num(16 + 7: end, 1));
    Vk_1 = car.axle_twist(g_gs_num_1, Jac_num_1, Jac_der_num_1, pb.x_1.*X_scale);  
    
    % cost function 
    fdtdak = 1/(pb.x(6)*X_scale(6) + tol);
    % Use this cost function to obtain Meccanica results
    %Jj = (fdtdak*dalfa)^2 + 10^-6*((Vk(2) - Vk_1(2))/dalfa)^2 + 1*(pb.u(7) - pb.u_1(7))^2 + 0.1*(pb.u(6) - pb.u_1(6))^2 + 10^-4*((Vk(6) - Vk_1(6))/dalfa)^2 + (pb.x_1(9)*X_scale(9))^2 + activation_comp*(pb.z(7)*pb.z(8))^2;
    Jj = (fdtdak*dalfa)^2 + 10^-2*(Vk(2) - Vk_1(2))^2 + 1*(pb.u(7) - pb.u_1(7))^2 + 0.1*(pb.u(6) - pb.u_1(6))^2 + 10^-4*(Vk(6) - Vk_1(6))^2 + ((pb.x_1(3) - pb.x(3))*X_scale(3))^2;   
    J_in = (Vk(1)/20 - 1)^2 + pb.x(1)^2 + 0.1*pb.x(6)^2;
    pb.Jj = activation_opt*Jj + activation_start*J_in;
    
    % comment this line to obtain Meccanica results
    pb.append_g(activation_comp*(pb.z(7)*pb.z(8)), -1e-4, 1e-4);
    
    % Stay on track
    pb.append_g(pb.x(1), -(pb.pgb(1)), pb.pgb(2));
    
    % Fzij constraints
    wrenchk = car.axle_wrench_sym(g_gs_num_k, Jac_num_k, Jac_der_num_k, pb.x.*X_scale, pb.u(1:end-1).*U_scale(1:end-1)); 
    Dz1_k  = DZ_1(wrenchk(4), pb.z(5)*Z_scale(5), pb.z(6)*Z_scale(6));
    Dz2_k  = DZ_2(wrenchk(4), pb.z(5)*Z_scale(5), pb.z(6)*Z_scale(6));
    Fzaero1_k = 0.5*rho*cz1*S*Vk(1)^2;
    Fzaero2_k = 0.5*rho*cz2*S*Vk(1)^2;
    Fz_11_k = Fz11tyre(wrenchk(3), Fzaero1_k, Fzaero2_k, wrenchk(5), Dz1_k); 
    Fz_12_k = Fz12tyre(wrenchk(3), Fzaero1_k, Fzaero2_k, wrenchk(5), Dz1_k); 
    Fz_21_k = Fz21tyre(wrenchk(3), Fzaero1_k, Fzaero2_k, wrenchk(5), Dz2_k);
    Fz_22_k = Fz22tyre(wrenchk(3), Fzaero1_k, Fzaero2_k, wrenchk(5), Dz2_k);
    pb.append_g((pb.z(1:4) - [Fz_11_k; Fz_12_k; Fz_21_k; Fz_22_k]./Z_scale(1:4)), zeros(4,1), zeros(4,1));
    
    % Constitutive constraint on Fx
    Fxa_k = pb.z(7)*Fx_scale;
    Fxb_k = pb.z(8)*Fx_scale;
    pb.append_g((pb.z(7) + pb.z(8)) - pb.u(1), -0., 0);

    % Adherence constraint 11
    Fy_11_k = Fytyre(alfa1_fun(Vk, pb.u.*U_scale), pb.z(1)*Z_scale(1));
    Fymax_11_k = Fymaxij(pb.z(1)*Z_scale(1));
    Fx_11_k = Fxb_k*kb/2;
    Fxmax_11_k = Fxmaxij(pb.z(1)*Z_scale(1));
    pb.append_g(((Fy_11_k/Fymax_11_k)^2 + (Fx_11_k/Fxmax_11_k)^2 - 1), -3, 0);
   
    % Adherence constraint 12
    Fy_12_k = Fytyre(alfa1_fun(Vk, pb.u.*U_scale), pb.z(2)*Z_scale(2));
    Fymax_12_k = Fymaxij(pb.z(2)*Z_scale(2));
    Fx_12_k = Fxb_k*kb/2;
    Fxmax_12_k = Fxmaxij(pb.z(2)*Z_scale(2));
    pb.append_g(((Fy_12_k/Fymax_12_k)^2 + (Fx_12_k/Fxmax_12_k)^2 - 1), -3, 0);
        
    % Adherence constraint 21
    Fy_21_k = Fytyre(alfa2_fun(Vk, pb.u.*U_scale), pb.z(3)*Z_scale(3));
    Fymax_21_k = Fymaxij(pb.z(3)*Z_scale(3));
    Fx_21_k = Fxa_k/2 + Fxb_k*(1-kb)/2;
    Fxmax_21_k = Fxmaxij(pb.z(3)*Z_scale(3));
    pb.append_g(((Fy_21_k/Fymax_21_k)^2 + (Fx_21_k/Fxmax_21_k)^2 - 1), -3, 0);
        
    % Adherence constraint 22
    Fy_22_k = Fytyre(alfa2_fun(Vk, pb.u.*U_scale), pb.z(4)*Z_scale(4));
    Fymax_22_k = Fymaxij(pb.z(4)*Z_scale(4));
    Fx_22_k = Fxa_k/2 + Fxb_k*(1-kb)/2;
    Fxmax_22_k = Fxmaxij(pb.z(4)*Z_scale(4));
    pb.append_g(((Fy_22_k/Fymax_22_k)^2 + (Fx_22_k/Fxmax_22_k)^2 - 1), -3, 0);
    
    % Constitutive constraints on Fy & Mz
    Fy_tot_k = Fy_11_k + Fy_12_k + Fy_21_k + Fy_22_k;
    Mz_tot_k = (Fy_11_k + Fy_12_k)*car.a1 - (Fy_21_k + Fy_22_k)*car.a2;
    pb.append_g(([pb.u(2); pb.u(6)] - [Fy_tot_k; Mz_tot_k]./[Fy_scale; Mz_scale]), zeros(2,1), zeros(2,1));    
    
    % Yi constraints
    pb.append_g((pb.z(5:6) - [Fy_11_k + Fy_12_k; Fy_21_k + Fy_22_k]./Z_scale(5:6)), zeros(2,1), zeros(2,1));
    
    % Power Limit
    pb.append_g(Vk(1)*pb.u(1)*Fx_scale/abs(Pmin),Pmin/abs(Pmin), Pmax/abs(Pmin));
   
    
    % Build map function 
    pb.build_map('JN', N-1);
    % Build g 
    pb.build_g;
    % Build J (different cost function for initial step)
    pb.build_J('range', 2:N);
    
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
    
    %%% add cost function term for initial step %%%
    kd = 0.001;
    fdtdak = 1/(pb.X(6,1)*X_scale(6) + tol);
    V1 = car.axle_twist(g_gs_num(:, :, 2), Jac_num(:, 2), Jac_der_num(:, 2), pb.X(:,1).*X_scale);
    V0 = car.axle_twist(g_gs_num(:, :, 1), Jac_num(:, 1), Jac_der_num(:, 1), pb.X_1(:,1).*X_scale);
    dalfa = dalfa_vec(1);
    pb.J = pb.J + activation_opt*((fdtdak*dalfa)^2 + 10*pb.X_1(10,1)^2 + 10*pb.X_1(11,1)^2 + 10*pb.X(10,1)^2 + 10*pb.X(11,1)^2 + 100*kd*(V0(2)^2 + V0(6)^2 + 100*pb.X_1(7,1)^2) + (pb.X_1(5,1)*X_scale(5))^2 + 10*(pb.X_1(4,1)*X_scale(4))^2 + 10*(pb.X_1(3,1)*X_scale(3))^2 + (pb.X(9,1)*X_scale(9))^2 + (pb.X_1(9,1)*X_scale(9))^2 + 10^-7*((V1(2) - V0(2))/dalfa)^2 + 10^-4*((V1(6) - V0(6))/dalfa)^2);
    pb.J = pb.J + activation_start*((V0(1)/20 - 1)^2 + (V1(1)/20 - 1)^2 + pb.X(1,1)^2 + 0.1*pb.X(6,1)^2);
    
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
        Zcons_h = casadi.(sym_type).sym(['Zcons' num2str(k_head)], o*(nx + nu + nz) + nx);
        Zcons = [Zcons; Zcons_h];        
        Ycons_h = casadi.(sym_type).sym(['Ycons' num2str(k_head)], o*(nx + nu + nz) + nx);
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
        Zcons_t = casadi.(sym_type).sym(['Zcons' num2str(k_tail)], o*(nx + nu + nz) + nx);
        Zcons = [Zcons; Zcons_t];
        Ycons_t = casadi.(sym_type).sym(['Ycons' num2str(k_tail)], o*(nx + nu + nz) + nx);
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
    
    
    problema = solver;
    numericalData.w0 = pb.w0;
    numericalData.initialCons = zeros(length(Zcons), 1);
    numericalData.lbw = pb.lbw;
    numericalData.ubw = pb.ubw;
    numericalData.lbg = pb.lbg;
    numericalData.ubg = pb.ubg;
    numericalData.nx = nx;
    numericalData.nu = nu;
    numericalData.nz = nz;
    numericalData.opts = opts;
    numericalData.d = d;
    
end