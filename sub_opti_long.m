function [problema, numericalData, scale]  = sub_opti_long(alfarange, pista, overlap_tail, overlap_head, problem_number, direct_c, init_subrange, opts)
    
% in this version we overlap only one set of states at the head and the
% tail (no controls, no algebraic variables)

    %% Discretization
    global alfa_end
    dalfa = round(alfarange(2) - alfarange(1),15);
    N = length(alfarange)-1;
    activation = MX.sym('act');
    activation_comp = MX.sym('act_comp');
    activation_opt = MX.sym('act_opt');
    activation_start = MX.sym('act_start');
    RHO_head = MX.sym('RHO_head');
    RHO_tail = MX.sym('RHO_tail');
    
    %% Other script
    import casadi.*
    if direct_c 
        direct_collocation;
        tau_colloc = tau_root(2:end);
    end
    car_parameters_ocp;
    common_ocp;
%     X_lb(1) = max(0,alfarange(1)-10*dalfa);
%     X_ub(1) = min(alfa_end, alfarange(end)+10*dalfa);
    
    %% Import initial guess
    guess_0 = init_subrange;
    guess_0.q = guess_0.q(:, 2:end);

    %% Construct NLP

    % Start with an empty NLP
                      
    J  = 0;
    Jc = 0;
    J_in = 0;
    prealloc_dim = floor(N*(nx + nu + nz)*1.5);
    w    = CellPrealloc(prealloc_dim);
    w0   = CellPrealloc(prealloc_dim);
    lbw  = CellPrealloc(prealloc_dim);
    ubw  = CellPrealloc(prealloc_dim);
    g    = CellPrealloc(prealloc_dim);
    lbg  = CellPrealloc(prealloc_dim);
    ubg  = CellPrealloc(prealloc_dim);
    Zw = CellPrealloc(prealloc_dim);
    Yw = CellPrealloc(prealloc_dim);
    Z0 = CellPrealloc(prealloc_dim); % initial consensus vector
    
    %% Definition of initial conditions (alfa = alfa_0)
    X_0 = guess_0.q(1,:)'./X_scale; 
    X0  = MX.sym('X0', nx);
    w.append(X0);
    
    if problem_number > 1
        %%%%% In opti singolo ef = 0 %%%%%%%
        lbw.append(X_lb);
        ubw.append(X_ub);
    else
        lbw.append([X_lb(1); 0; X_lb(3:end)]);
        ubw.append([X_ub(1); 0; X_ub(3:end)]);
    end
    
    w0.append(X_0);
    Xk = X0;
    
    %Initial constraint
    Vk = car.axle_twist(g_gs_num(:, :, 1), Jac_num(:, 1), Jac_der_num(:, 1), Xk.*X_scale); 
    if problem_number > 1
        g.append(Vk(1)/Vmax);
        lbg.append(0);
        ubg.append(1);%((Pmax/(0.5*rho*S*cx))^(1/3))
    else
        g.append(Vk(1)/Vmax);
        lbg.append(Vi/Vmax);
        ubg.append(Vi/Vmax);
    end
        
    %%%%%%  %%%%%%
    pos = pista.fun_pos(alfarange(1));
    g.append(Xk(1));
    lbg.append(-(pos(5)-t/2)/ep_scale);
    ubg.append( (pos(5)-t/2)/ep_scale);
    

    guess_0.q(:,1) = 0;
    guess_0.q(:,7) = 0;
    
    %% Loop to build NLP
    for k = 1:N

        %Define previous variables
        if k == 1
            Xk_1 = Xk;
        else
            Xk_1 = Xk;
            Zk_1 = Zk;
            Uk_1 = Uk;
        end
                
        %Initialization of states
        X_0 = guess_0.q(k+1,:)'./X_scale;%scaled_states_0;
        %Initialization of controls
        U_0 = guess_0.w(k,:)'./U_scale(1:nu-1);%scaled_controls_0;
        % Initialization for steer angle (Fy11/Cf + beta)
        Vk_0 = car.axle_twist(g_gs_num(:, :, k+1), Jac_num(:, k+1), Jac_der_num(:, k+1), X_0.*X_scale);
        alfa1num = (U_0(2)*U_scale(2)*a2 + U_0(6)*U_scale(6))/(2*(a1+a2))/CsFz(m*G*a2/(a1+a2)/2);
        delta0 = alfa1num + (Vk_0(2) + Vk_0(end)*a1)/Vk_0(1);
        U_0 = [U_0;delta0/U_scale(nu)];
        
        %Initialization of algebraic parameters
        wrenchout_0 = car.axle_wrench_sym(g_gs_num(:, :, k+1), Jac_num(:, k+1), Jac_der_num(:, k+1),X_0.*X_scale, U_0(1:end-1).*U_scale(1:end-1));
        Dz1_0  = DZ_1(wrenchout_0(4), U_0(2)*U_scale(2)*(a2/l) + U_0(6)*U_scale(6)/l, U_0(2)*U_scale(2)*(a1/l) - U_0(6)*U_scale(6)/l);
        Dz2_0  = DZ_2(wrenchout_0(4), U_0(2)*U_scale(2)*(a2/l) + U_0(6)*U_scale(6)/l, U_0(2)*U_scale(2)*(a1/l) - U_0(6)*U_scale(6)/l);
        Fzaero1_0 = 0.5*rho*cz1*S*Vk_0(1)^2;
        Fzaero2_0 = 0.5*rho*cz2*S*Vk_0(1)^2;
        Fz11_0 = Fz11tyre(wrenchout_0(3), Fzaero1_0, Fzaero2_0, wrenchout_0(5), Dz1_0);
        Fz12_0 = Fz12tyre(wrenchout_0(3), Fzaero1_0, Fzaero2_0, wrenchout_0(5), Dz1_0);
        Fz21_0 = Fz21tyre(wrenchout_0(3), Fzaero1_0, Fzaero2_0, wrenchout_0(5), Dz2_0);
        Fz22_0 = Fz22tyre(wrenchout_0(3), Fzaero1_0, Fzaero2_0, wrenchout_0(5), Dz2_0);
        Y1_0 = Fytyre(alfa1_fun(Vk_0, U_0.*U_scale), Fz11_0)+Fytyre(alfa1_fun(Vk_0, U_0.*U_scale), Fz12_0);
        Y2_0 = Fytyre(alfa2_fun(Vk_0, U_0.*U_scale), Fz21_0)+Fytyre(alfa2_fun(Vk_0, U_0.*U_scale), Fz22_0);
        Fxa_0 = max(U_0(1)*U_scale(1), 0);
        Fxb_0 = min(U_0(1)*U_scale(1), 0);
        Z_0 = [Fz11_0; Fz12_0; Fz21_0; Fz22_0; Y1_0; Y2_0; Fxa_0; Fxb_0]./Z_scale;%scaled_algebraic_0;
        
        %Consensus for first iteration
        if k == floor(overlap_head/2)  && overlap_head > 0 
            Z0.append(X_0);
        elseif k == 1 && floor(overlap_head/2) + 1 == 1 && overlap_head > 0 
            Z0.append(X_0);
        end
        if k == floor(overlap_head/2) + 1  && overlap_head > 0 
            Z0.append([U_0; Z_0; X_0]);
        end
        
        
        if k== N -(floor(overlap_tail/2)+1) && overlap_tail > 0
            Z0.append(X_0);
        end
        
        if k== N-floor(overlap_tail/2) && overlap_tail > 0
            Z0.append([U_0; Z_0; X_0]);
        end
        
        % Define k-th controls vector
        Uk = MX.sym(['U_' num2str(k)], nu);
        w.append(Uk);
        lbw.append(U_lb);
        ubw.append(U_ub);
        w0.append(U_0);


        % Define k-th algebraic vector
        Zk = MX.sym(['Z_' num2str(k)], nz);
        w.append(Zk);
        lbw.append(Z_lb);
        ubw.append(Z_ub);
        w0.append(Z_0);
    
        if direct_c
            % State at collocation points
            Xkj = cell(d, 1);
            for j=1:d
                Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], nx);
                w.append(Xkj{j});
                lbw.append(X_lb);
                ubw.append(X_ub);
                w0.append(X_0);
            end

            % Loop over collocation points
            Xk_end = D(1)*Xk;
            for j=1:d
               % Expression for the state derivative at the collocation point
               xp = C(1,j+1)*Xk;
               for r=1:d
                   xp = xp + C(r+1,j+1)*Xkj{r};
               end

               % Append collocation equations
               g_gs_colloc = full(pista.fun_g_gs(alfarange(k) + dalfa*tau_colloc(j)));
               Jac_colloc = full(car.T_jac(alfarange(k) + dalfa*tau_colloc(j)));
               Jac_der_colloc = full(car.T_jac_der(alfarange(k) + dalfa*tau_colloc(j)));
               fj = F(g_gs_colloc, Jac_colloc, Jac_der_colloc, Xkj{j}.*X_scale, Uk.*U_scale);
               g.append((dalfa*fj(2:end) - xp.*X_scale)./X_scale);
               lbg.append(zeros(nx,1));
               ubg.append(zeros(nx,1));

               % Add contribution to the end state
               Xk_end = Xk_end + D(j+1)*Xkj{j};

               % Add contribution to quadrature function
               % J = J + B(j+1)*qj*dalfa;
            end
        end

        % Define (k+1)-th states vector
        Xk = MX.sym(['X_' num2str(k+1)], nx);
        w.append(Xk);
        lbw.append(X_lb);
        ubw.append(X_ub);
        w0.append(X_0);
        
        if direct_c
            g.append((Xk_end(1:end)-Xk(1:end)));
            lbg.append(zeros(nx,1));
            ubg.append(zeros(nx,1));  
        end
        
        if ~direct_c
            % Continuity constraints (Integration scheme)
            g.append(Xk_dot*dalfa./X_scale - (Xk - Xk_1));
            lbg.append(zeros(node,1));
            ubg.append(zeros(node,1));
        end

        
        % Previous twist
        Vk_1 = Vk;

        % Actual twist
        Vk = car.axle_twist(g_gs_num(:, :, k+1), Jac_num(:, k+1), Jac_der_num(:, k+1), Xk.*X_scale);  
        
        % Evaluate derivative at collocation point
        Xk_dot = F(g_gs_num(:, :, k+1), Jac_num(:, k+1), Jac_der_num(:, k+1), Xk.*X_scale,  Uk.*U_scale);

        % Evaluate cost function at k-step
        fdtdak = 1/(Xk(6)*X_scale(6)+tol); 
        kd = 0.001;
        
        if k == 1
            J = J + (fdtdak*dalfa)^2 + Xk_1(10)^2 + Xk_1(11)^2 + Xk(10)^2 + Xk(11)^2 + 100*kd*(Vk_1(2)^2 + Vk_1(6)^2 + 100*Xk_1(7)^2) + (Xk_1(9)*X_scale(9))^2 + 10^-7*((Vk(2) - Vk_1(2))/dalfa)^2 + 10^-4*((Vk(6) - Vk_1(6))/dalfa)^2;
            J_in = J_in + (Vk_1(1)/20 - 1)^2 + (Vk(1)/20 - 1)^2 + Xk(1)^2 + 0.1*Xk(6)^2;
            if floor(overlap_head/2) + 1 == 1 && overlap_head > 0
                
                % consensus term
                Zcons = MX.sym(['Zcons' num2str(k)], 2*nx + nu + nz);
                Zw.append(Zcons);
                Ycons = MX.sym(['Ycons' num2str(k)], 2*nx + nu + nz);
                Yw.append(Ycons);
                cdiff = ([Xk_1; Uk; Zk; Xk] - Zcons);
                Jc = Jc + Ycons'*cdiff + 0.5*RHO_head*(cdiff.')*(cdiff);
            end
        else
            J = J + (fdtdak*dalfa)^2 + 10^-2*((Vk(2) - Vk_1(2)))^2 + 1*(Uk(7)-Uk_1(7))^2 + 0.1*(Uk(6)-Uk_1(6))^2 + 10^-4*((Vk(6) - Vk_1(6)))^2 + ((Xk_1(3) - Xk_1(3))*X_scale(9))^2;
            J_in = J_in + (Vk(1)/20 - 1)^2 + Xk(1)^2 +0.1*Xk(6)^2;
            
            if k == floor(overlap_head/2) + 1 && overlap_head > 0 % head consensus
                % consensus term
                Zcons = MX.sym(['Zcons' num2str(k)], 2*nx + nu + nz);
                Zw.append(Zcons);
                Ycons = MX.sym(['Ycons' num2str(k)], 2*nx + nu + nz);
                Yw.append(Ycons);
                cdiff = ([Xk_1; Uk; Zk; Xk] - Zcons);
                Jc = Jc + Ycons'*cdiff + 0.5*RHO_head*(cdiff.')*(cdiff);
            end
            
            if k == N-floor(overlap_tail/2) && overlap_tail > 0 % tail consensus
                % consensus term
                Zcons = MX.sym(['Zcons' num2str(k)], 2*nx + nu + nz);
                Zw.append(Zcons);
                Ycons = MX.sym(['Ycons' num2str(k)], 2*nx + nu + nz);
                Yw.append(Ycons);
                cdiff = ([Xk_1; Uk; Zk; Xk] - Zcons);
                Jc = Jc + Ycons'*cdiff + 0.5*RHO_tail*(cdiff.')*(cdiff);
            end
            
        end

    % Stay on track constraint
    pos = pista.fun_pos(alfarange(k+1));
    g.append(Xk(1));
    lbg.append( -(pos(5)-t/2)/ep_scale);
    ubg.append((pos(5)-t/2)/ep_scale);
        
    g.append(activation_comp*(Zk(7)*Zk(8))^2);
    lbg.append(-1e-6);
    ubg.append(1e-6);
    % Evaluate output wrench
    wrenchk = car.axle_wrench_sym(g_gs_num(:, :, k+1), Jac_num(:, k+1), Jac_der_num(:, k+1), Xk.*X_scale, Uk(1:end-1).*U_scale(1:end-1)); 
    Dz1_k  = DZ_1(wrenchk(4), Zk(5)*Z_scale(5), Zk(6)*Z_scale(6));
    Dz2_k  = DZ_2(wrenchk(4), Zk(5)*Z_scale(5), Zk(6)*Z_scale(6));
    Fzaero1_k = 0.5*rho*cz1*S*Vk(1)^2;
    Fzaero2_k = 0.5*rho*cz2*S*Vk(1)^2;
    Fz_11_k = Fz11tyre(wrenchk(3), Fzaero1_k, Fzaero2_k, wrenchk(5), Dz1_k); 
    Fz_12_k = Fz12tyre(wrenchk(3), Fzaero1_k, Fzaero2_k, wrenchk(5), Dz1_k); 
    Fz_21_k = Fz21tyre(wrenchk(3), Fzaero1_k, Fzaero2_k, wrenchk(5), Dz2_k);
    Fz_22_k = Fz22tyre(wrenchk(3), Fzaero1_k, Fzaero2_k, wrenchk(5), Dz2_k);
    % Fzij constraints
    g.append((Zk(1:4) - [Fz_11_k; Fz_12_k; Fz_21_k; Fz_22_k]./Z_scale(1:4)));
    lbg.append(zeros(4,1));
    ubg.append(zeros(4,1));
    
    
    Fxa_k = Zk(7)*Fx_scale;
    Fxb_k = Zk(8)*Fx_scale;
    g.append((Zk(7) + Zk(8)) - Uk(1));
    lbg.append(-0.);
    ubg.append( 0.);    
    
    % Adherence constraint 11
    Fy_11_k = Fytyre(alfa1_fun(Vk, Uk.*U_scale), Zk(1)*Z_scale(1));
    Fymax_11_k = Fymaxij(Zk(1)*Z_scale(1));
    Fx_11_k = Fxb_k*kb/2;
    Fxmax_11_k = Fxmaxij(Zk(1)*Z_scale(1));
    g.append(((Fy_11_k/Fymax_11_k)^2 + (Fx_11_k/Fxmax_11_k)^2 - 1));
    lbg.append( -3.);
    ubg.append(0.);
    
    % Adherence constraint 12
    Fy_12_k = Fytyre(alfa1_fun(Vk, Uk.*U_scale), Zk(2)*Z_scale(2));
    Fymax_12_k = Fymaxij(Zk(2)*Z_scale(2));
    Fx_12_k = Fxb_k*kb/2;
    Fxmax_12_k = Fxmaxij(Zk(2)*Z_scale(2));
    g.append(((Fy_12_k/Fymax_12_k)^2 + (Fx_12_k/Fxmax_12_k)^2 - 1));
    lbg.append( -3.);
    ubg.append(0.);
    
    % Adherence constraint 21
    Fy_21_k = Fytyre(alfa2_fun(Vk, Uk.*U_scale), Zk(3)*Z_scale(3));
    Fymax_21_k = Fymaxij(Zk(3)*Z_scale(3));
    Fx_21_k = Fxa_k/2 + Fxb_k*(1-kb)/2;
    Fxmax_21_k = Fxmaxij(Zk(3)*Z_scale(3));
    g.append(((Fy_21_k/Fymax_21_k)^2 + (Fx_21_k/Fxmax_21_k)^2 - 1));
    lbg.append(-3.);
    ubg.append(0.);
    
    % Adherence constraint 22
    Fy_22_k = Fytyre(alfa2_fun(Vk, Uk.*U_scale), Zk(4)*Z_scale(4));
    Fymax_22_k = Fymaxij(Zk(4)*Z_scale(4));
    Fx_22_k = Fxa_k/2 + Fxb_k*(1-kb)/2;
    Fxmax_22_k = Fxmaxij(Zk(4)*Z_scale(4));
    g.append(((Fy_22_k/Fymax_22_k)^2 + (Fx_22_k/Fxmax_22_k)^2 - 1));
    lbg.append(-3.);
    ubg.append(0.); 
        
    % Constitutive constraints
    Fy_tot_k = Fy_11_k + Fy_12_k + Fy_21_k + Fy_22_k;
    Mz_tot_k = (Fy_11_k + Fy_12_k)*a1 - (Fy_21_k + Fy_22_k)*a2;
    g.append(([Uk(2); Uk(6)] - [Fy_tot_k; Mz_tot_k]./[Fy_scale; Mz_scale]));
    lbg.append(zeros(2,1));
    ubg.append(zeros(2,1));
    
    
    % Yi constraints
    g.append((Zk(5:6) - [Fy_11_k + Fy_12_k; Fy_21_k + Fy_22_k]./Z_scale(5:6)));
    lbg.append(zeros(2,1));
    ubg.append(zeros(2,1));

    
    % Power Limit
    g.append(Vk(1)*Uk(1)*Fx_scale/abs(Pmin));
    lbg.append(Pmin/abs(Pmin));
    ubg.append(Pmax/abs(Pmin));

    end
        
%     g = [g(:)', {Xk(3)}];
%     lbg = [lbg; -0.01/X_scale(3)];
%     ubg = [ubg;  0.01/X_scale(3)];

    J = activation_opt*J + activation.*(Jc) + activation_start*J_in;
    
    % Create an NLP solver 
    prob = struct('f', J, 'x', w.flattenData, 'g', g.flattenData, 'p', vertcat(Zw.flattenData, Yw.flattenData, activation, RHO_head, RHO_tail, activation_comp, activation_opt, activation_start));
    
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
    numericalData.w0 = w0.flattenData;
    numericalData.initialCons = Z0.flattenData;
    numericalData.lbw = lbw.flattenData;
    numericalData.ubw = ubw.flattenData;
    numericalData.lbg = lbg.flattenData;
    numericalData.ubg = ubg.flattenData;
    numericalData.nx = nx;
    numericalData.nu = nu;
    numericalData.nz = nz;
    numericalData.opts = opts;
%     numericalData.dll_filename = sprintf('%s_%i.dll', 'Compilati\sub_problem', problem_number);
%     numericalData.bat_filename = sprintf('%s_%i.bat', 'Compilati\compilation', problem_number);  
    if direct_c 
        numericalData.d = d;
    end
    
end