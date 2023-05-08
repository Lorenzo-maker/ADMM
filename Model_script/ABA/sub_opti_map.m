function [problema, numericalData, scale]  = sub_opti_map(alfarange, pista, o, id, ID, problem_number, d, init_subrange, alpha_vec, opts)
    
% in this version we overlap only one set of states at the head and the
% tail (no controls, no algebraic variables)

    %% Discretization
    dalfa = round(alfarange(2) - alfarange(1),15);
    N = length(alfarange)-1;
    
    %% Other script
    import casadi.*
    direct_collocation;    
    car_parameters_ocp;
    common_ocp;

    %% Import initial guess
    guess_0 = init_subrange;
    guess_0.q = guess_0.q(:, 2:end); % remove alpha from guess
    X_0 = guess_0.q(1,:)'./X_scale; 
    guess_0.q = guess_0.q;
%     X_0 = zeros(nx,1);
%     X_0(3) = -0.022;
%     X_0(6) = 6e-4;
%     X_0 = X_0./X_scale;
    
    %% Construct NLP
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
    np = 6;
    Npg = (4*4 + 6 + 6)*(d+2);
    Npgb = 2;
    Pdata = struct;
    Pdata.Pg = zeros((16 + 6 + 6)*(d+2), N);
    for k = 1:N
        Pdata.Pg(1:(16+6+6),k) = [reshape(g_gs_num(:, :, k),16,1); Jac_num(:, k); Jac_der_num(:, k)];
        for j = 1:d
            Pdata.Pg((16+6+6) + (j-1)*(16+6+6) + 1:(16+6+6) + j*(16+6+6),k) = [g_gs_colloc{k, j}(:); Jac_colloc{k, j}; Jac_der_colloc{k, j}];
        end
        Pdata.Pg(end-(16+6+6)+1:end,k) = [reshape(g_gs_num(:, :, k+1),16,1); Jac_num(:, k+1); Jac_der_num(:, k+1)];
    end
    pos = pista.fun_pos(alfarange);
    Pdata.Pgb = [(pos(5,2:end) - t/2); pos(5, 2:end) - t/2]./ep_scale;
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
    
%     guess_0.q(:,1) = 0;
%     guess_0.q(:,7) = 0;
    
    % Initial guess
    for k = 1:N
        
        V0 = Vi;
        X_0 = zeros(11,1);
        X_0(3) = mean(guess_0.q(:,3))/X_scale(3);%-0.022/X_scale(3);%
        X_0(6) = mean(guess_0.q(:,6))/X_scale(6);%6e-04/X_scale(6);%mean(guess_0.q(:,6))/X_scale(6);%dalfa/(ds/(V0))/X_scale(6);
        Fz0 = car.M(1,1)*9.81 + 0.5*rho*(cz1+cz2)*S*V0^2;
        U_0 = [0.5*rho*cx*S*V0^2; 0; car.M(1,1)*9.81 + 0.5*rho*(cz1+cz2)*S*V0^2; 0; 0; 0; 0]./U_scale;
        Z_0 = [Fz0/4*ones(4,1); 0; 0; 0.5*rho*cx*S*V0^2; 0]./Z_scale;
        pb.w0 = [pb.w0; U_0; Z_0; repmat(X_0, d+1, 1)];
    end
    pb.w0 = [guess_0.q(1,:)'./X_scale; pb.w0];%[X_0./X_scale; pb.w0]; %
    
    % Start constraint
    Vk = car.axle_twist(g_gs_num(:, :, 1), Jac_num(:, 1), Jac_der_num(:, 1), pb.x_1.*X_scale); 
    if problem_number > 1
        pb.append_g(Vk(1)/Vmax, 0, 1, 'start');
    else
        pb.append_g(Vk(1)/Vmax, Vi/Vmax, Vi/Vmax, 'start');               
    end
    pb.append_g(pb.x_1(1), -(pos(5,1)-t/2)/ep_scale, (pos(5,1)-t/2)/ep_scale, 'start');
     
    % Collocation
    data_num = reshape(pb.pg, 16+6+6, d+2);
    for i = 1:d
        g_gs_colloc_j = reshape(data_num(1:16, i+1),4,4);
        Jac_colloc_j = (data_num(16 + 1: 16 + 6, i+1));
        Jac_der_colloc_j = (data_num(16 + 7: end, i+1));
        fj = F(g_gs_colloc_j, Jac_colloc_j, Jac_der_colloc_j, pb.xc(:,i).*X_scale, pb.u.*U_scale);
        pb.append_g((dalfa*fj(2:end) - pb.xp(:,i).*X_scale)./X_scale, zeros(nx,1), zeros(nx, 1));
    end
    
    pb.append_g(pb.xc_end - pb.x, zeros(nx,1), zeros(nx, 1));
    
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
    %%%%%%%
    pb.append_g(pb.x(1), -(pb.pgb(1)), pb.pgb(2));
    
    wrenchk = car.axle_wrench_sym(g_gs_num_k, Jac_num_k, Jac_der_num_k, pb.x.*X_scale, pb.u(1:end-1).*U_scale(1:end-1)); 
    Dz1_k  = DZ_1(wrenchk(4), pb.z(5)*Z_scale(5), pb.z(6)*Z_scale(6));
    Dz2_k  = DZ_2(wrenchk(4), pb.z(5)*Z_scale(5), pb.z(6)*Z_scale(6));
    Fzaero1_k = 0.5*rho*cz1*S*Vk(1)^2;
    Fzaero2_k = 0.5*rho*cz2*S*Vk(1)^2;
    Fz_11_k = Fz11tyre(wrenchk(3), Fzaero1_k, Fzaero2_k, wrenchk(5), Dz1_k); 
    Fz_12_k = Fz12tyre(wrenchk(3), Fzaero1_k, Fzaero2_k, wrenchk(5), Dz1_k); 
    Fz_21_k = Fz21tyre(wrenchk(3), Fzaero1_k, Fzaero2_k, wrenchk(5), Dz2_k);
    Fz_22_k = Fz22tyre(wrenchk(3), Fzaero1_k, Fzaero2_k, wrenchk(5), Dz2_k);
    % Fzij constraints
    pb.append_g((pb.z(1:4) - [Fz_11_k; Fz_12_k; Fz_21_k; Fz_22_k]./Z_scale(1:4)), zeros(4,1), zeros(4,1));
    %
    Fxa_k = pb.z(7)*Fx_scale;
    Fxb_k = pb.z(8)*Fx_scale;
    pb.append_g((pb.z(7) + pb.z(8)) - pb.u(1), -0., 0);
    %%%%
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
    
    % Constitutive constraints
    Fy_tot_k = Fy_11_k + Fy_12_k + Fy_21_k + Fy_22_k;
    Mz_tot_k = (Fy_11_k + Fy_12_k)*a1 - (Fy_21_k + Fy_22_k)*a2;
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
    pb.J = pb.J + activation_opt*((fdtdak*dalfa)^2 + 10*pb.X_1(10,1)^2 + 10*pb.X_1(11,1)^2 + 10*pb.X(10,1)^2 + 10*pb.X(11,1)^2 + 100*kd*(V0(2)^2 + V0(6)^2 + 100*pb.X_1(7,1)^2) + (pb.X_1(5,1)*X_scale(5))^2 + 10*(pb.X_1(4,1)*X_scale(4))^2 + 10*(pb.X_1(3,1)*X_scale(3))^2 + (pb.X(9,1)*X_scale(9))^2 + (pb.X_1(9,1)*X_scale(9))^2 + 10^-7*((V1(2) - V0(2))/dalfa)^2 + 10^-4*((V1(6) - V0(6))/dalfa)^2);
    pb.J = pb.J + activation_start*((V0(1)/20 - 1)^2 + (V1(1)/20 - 1)^2 + pb.X(1,1)^2 + 0.1*pb.X(6,1)^2);
    
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
    numericalData.nx = nx;
    numericalData.nu = nu;
    numericalData.nz = nz;
    numericalData.opts = opts;
%     numericalData.dll_filename = sprintf('%s_%i.dll', 'Compilati\sub_problem', problem_number);
%     numericalData.bat_filename = sprintf('%s_%i.bat', 'Compilati\compilation', problem_number);  
    numericalData.d = d;
    
end