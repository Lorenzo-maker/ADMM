classdef nlp_casadi_2 < handle
    properties
        nx        % State dimension
        nu        % Control dimension
        nz        % Algebraic parameters dimension  
        nuc       % Control dimension in colloc points
        nzc       % Algebraic parameters dimension in colloc points
        nwg       % Global optimization variables          
        wg        % Symbolic vector of global variable
        x0        % Symbolic vector states in starting mesh node 
        x         % Symbolic vector states in k+1 mesh node
        x_1       % Symbolic vector states in k mesh node
        xc        % Symbolic matrix states in colloc points (k / k+1) (nx, d)
        xp        % Symbolic matrix of states derivative in colloc points (k / k+1) (nx, d)
        xc_end    % Symbolic terminal value of the states at k+1
        x_end     % Symbolic vector states in final mesh node
        u0        % Symbolic controls in starting mesh node
        u         % Symbolic controls in k interval
        u_1       % Symbolic controls in k-1 interval
        uc        % Symbolic controls in colloc points (k / k+1) (not necessary with the same dimension of nu)
        u_end     % Symbolic controls in final mesh node       
        z0        % Symbolic ALg. parameters in starting mesh node        
        z         % Symbolic Alg. parameters in k interval
        z_1       % Symbolic Alg. parameters in k-1 interval
        zc        % Symbolic Alg. parameters in colloc points (k / k+1) (not necessary with the same dimension of nz)
        z_end     % Symbolic Alg. parameters in final mesh node
        N         % Number of discretization interval
        d         % Degree for collocation (hence, number of collocation points)
        D         % Matrix for xc_end evaluation
        C         % Matrix for xp evaluation 
        B         % Matrix for quadrature evaluation   
        Wg        % Vector of global variable (nwg, 1)
        X         % Matrix of the states (nx, N)
        X_1       % Matrix of the states shifted of one element backward (nx, N)
        Xc        % Matrix of the states in colloc points (nx*d, N)
        U         % Matrix of the controls (nu, N)
        U_1       % Matrix of the controls shifted one element backward (nu, N) 
        Uc        % Matrix of the controls in colloc points 
        Z         % Matrix of Alg. parameters (nz, N)       
        Z_1       % Matrix of Alg. parameters shifted one element backward (nz, N)
        Zc        % Matrix of Alg. parameters in colloc points 
        w         % Vector of optimization variables 
        lbw       % Vector of lower bound (X0, U, Z, Xc, Uc, Zc, X)
        ubw       % Vector of upper bound (X0, U, Z, Xc, Uc, Zc, X)
        Flbw      % Function to evaluate parametric lowr bound on variables
        Fubw      % Function to evaluate parametric upper bound on variables
        w0        % Intial guess vector    
        ng        % Length of loop constraints
        gv        % Vector of loop constraints, written one time
        gv_s      % Vector of start constraints 
        lbg_s     % Lower bound for start constraints
        ubg_s     % Upper bound for start constraints
        gv_e      % Vector of end constraints
        lbg_e     % Lower bounf for end constraints
        ubg_e     % Upper bound for end constraints    
        gv_c      % Vector of periodic constraints (start-end)
        lbg_c     % Lower bound of periodic constraints
        ubg_c     % Upper bound of periodic constraints 
        gv_m      % Manual constraints
        lbg_m     % Lower bound on manual constraint
        ubg_m     % Upper bound on manual constraint
        g         % Total constraints vector 
        lbg       % Vector of constraints lower bound
        ubg       % Vector of constraints upper bound
        pg        % Symbolic parameters that enter in constraint evaluation (i.e. g(x, u, z, pg) = 0)
        Npg       % Number of parameters that enter in constraints
        Pg        % Numerical value of pg along the grid
        pgb       % Symbolic parameters that enter in lower or upper bound  (i.e. -pgb(1) <= g(x,u,z) <= pgb(2)) 
        Npgb      % Number of parametric lower bound + Number of parametric upper bound 
        Pgb       % Numerical value of pgb along the grid
        Eq_g      % Function to evaluate loop constraints 
        Eq_g_map  % Map version of Eq_g
        Eq_g_s    % Function to evaluate start constraints
        Eq_g_e    % Function to evaluate end constraints
        Eq_g_c    % Function to evaluate periodic constraints
        Eq_bg     % Function to evaluate lower and upper bound for loop constraints (it is function of some numerical parameters to have variable constraints in the grid)
        Eq_bg_map % Map version of Eq_bg       
        Jk        % Symbolic value of cost function in k interval
        J_fun     % Functio to evaluate Jk
        J_fun_map % Map version of J_fun
        Jj        % Vector of cost function evaluated at colloc points
        J         % Cost function
        Nthread   % Number of threads for Map function
        p         % Symbolic parameters generalized NLP
        np        % Number of symbolic parameters 
        sym_type  % Type for symbolic variables
        solver    % IPOPT solver
        sol       % NLP solution
        sol_num   % Numerical solution without first nx states
        X0_num    % Numerical solution for the first nx states
        w_global  % Numerical solution for global variables
        check_x0  % value of Xb.x0 
        colloc_type
    end
    
    methods
        function obj = nlp_casadi_2(nx, nu, nz, nuc, nzc, nwg, N, d, Xb, Ub, Zb, Ucb, Zcb, Wgb, Npg, Npgb, Pdata, Nthread, np, sym_type, colloc_type)
            arguments
                nx = 1; nu = 1; nz = 1; nuc = 0; nzc = 0; nwg = 0; N = 1; d = 2; Xb = []; Ub = []; Zb = []; ...
                Ucb = []; Zcb = []; Wgb = []; Npg = 0; Npgb = 0; Pdata = []; Nthread = 1; np = 0; sym_type = 'SX', colloc_type = 'legendre'; 
            end
            import casadi.*            
            obj.sym_type = sym_type;
            obj.colloc_type = colloc_type;
            obj.Nthread = Nthread;
            obj.N = N;
            obj.d = d;           
            % Adjust empty quantities
            Xb = obj.empty_struct(Xb,'x0','yes');
            Ub = obj.empty_struct(Ub);
            Zb = obj.empty_struct(Zb);
            Ucb = obj.empty_struct(Ucb);
            Zcb = obj.empty_struct(Zcb);
            Wgb = obj.empty_struct(Wgb);
            % Symbolical value on generic step grid
            obj.np = np;
            obj.nx = nx;
            obj.nu = nu;
            obj.nz = nz;
            obj.nuc = nuc;
            obj.nzc = nzc;
            obj.nwg = nwg;
            obj.Npg = Npg;
            obj.Npgb = Npgb;
            %
            obj.casadi_sym()
            X0 = obj.w(1:obj.nx);
            wr = reshape(obj.w(obj.nx+1:end-obj.nwg), (obj.nx + obj.nuc + obj.nzc)*d + obj.nu + obj.nz + obj.nx, obj.N);
            obj.U = wr(1:obj.nu, :);
            obj.U_1 = [obj.U(:,1), obj.U(:,1:end-1)];
            obj.Z = wr(obj.nu+1:obj.nu+obj.nz,:);
            obj.Z_1 = [obj.Z(:,1), obj.Z(:,1:end-1)];
            Xcc = wr(obj.nu+obj.nz+1:obj.nu+obj.nz+obj.nx*d, :);
            obj.Xc = reshape(Xcc, obj.nx, obj.N*d);
            Ucc = wr(obj.nu+obj.nz+obj.nx*d+1:obj.nu+obj.nz+(obj.nx+obj.nuc)*d, :);
            obj.Uc = reshape(Ucc, obj.nuc, obj.N*d);
            Zcc = wr(obj.nu+obj.nz + (obj.nx+obj.nuc)*d +1:obj.nu+obj.nz+(obj.nx+obj.nuc+obj.nzc)*d, :);
            obj.Zc = reshape(Zcc, obj.nzc, obj.N*d);
            obj.X = wr(end-obj.nx + 1: end,:);
            obj.X_1 = [X0, obj.X(:,1:end-1)];
            obj.Wg = obj.w(end - obj.nwg + 1:end);
            if ~isfield(Pdata,'Pg')
                obj.Pg = []; % Matrice con N colonne
            else
                obj.Pg = Pdata.Pg; % Matrice con N colonne
            end
            if ~isfield(Pdata,'Pgb')    
                obj.Pgb = []; % Matrice con N colonne con righe lower e upper
            else
                obj.Pgb = Pdata.Pgb; % Matrice con N colonne con righe lower e upper
            end
            % Other useful quantites
            obj.g = [];
            obj.col_points()
            obj.Jj = [];
            % Build lbw and ubw
            obj.lbw = [Ub.lb; Zb.lb; repmat(Xb.lb, d, 1); repmat(Ucb.lb, d, 1); repmat(Zcb.lb, d, 1); Xb.lb];
            obj.lbw = repmat(obj.lbw, obj.N, 1);
            obj.lbw = [obj.lbw; Wgb.lb];
            obj.ubw = [Ub.ub; Zb.ub; repmat(Xb.ub, d, 1); repmat(Ucb.ub, d, 1); repmat(Zcb.ub, d, 1); Xb.ub];
            obj.ubw = repmat(obj.ubw, obj.N, 1);
            obj.ubw = [obj.ubw; Wgb.ub];
            obj.check_x0 = Xb.x0;
            if ~isempty(Xb.x0) && isa(Xb.x0, 'double') %~isempty(Xb.x0)
                if ~isfield(Xb, 'x0_lb')
                    obj.lbw = [Xb.x0; obj.lbw];
                    obj.ubw = [Xb.x0; obj.ubw];
                else
                    obj.lbw = [Xb.x0_lb; obj.lbw];
                    obj.ubw = [Xb.x0_ub; obj.ubw];
                end
            elseif strcmp(Xb.x0, 'parametric')
                if strcmp(sym_type, 'SXMX')
                    x0p = SX.sym('x0p', obj.nx);
                else
                    x0p = casadi.(sym_type).sym('x0p', obj.nx);
                end
                obj.lbw = [x0p; obj.lbw];
                obj.ubw = [x0p; obj.ubw];
                obj.lbw = Function('Flbw', {x0p}, {full(obj.lbw)});
                obj.ubw = Function('Fubw', {x0p}, {full(obj.ubw)});
            else
                obj.lbw = [Xb.lb; obj.lbw];
                obj.ubw = [Xb.ub; obj.ubw];
            end           
            obj.w0 = [];
        end
        % Method to append constraint
        function append_g(obj, g, glb, gub, type)
            arguments
                obj;
                g;
                glb;
                gub;
                type = 'loop';
            end
            import casadi.*
            if strcmp(type, 'start')
                obj.gv_s = [obj.gv_s; g];
                obj.lbg_s = [obj.lbg_s; glb];
                obj.ubg_s = [obj.ubg_s; gub];
            elseif strcmp(type, 'end')
                obj.gv_e = [obj.gv_e; g];
                obj.lbg_e = [obj.lbg_e; glb];
                obj.ubg_e = [obj.ubg_e; gub];
            elseif strcmp(type, 'cyclic')
                obj.gv_c = [obj.gv_c; g];
                obj.lbg_c = [obj.lbg_c; glb];
                obj.ubg_c = [obj.ubg_c; gub];
            elseif strcmp(type, 'manually')
                obj.gv_m = [obj.gv_m; g];
                obj.lbg_m = [obj.lbg_m; glb];
                obj.ubg_m = [obj.ubg_m; gub];                
            else
                obj.gv = [obj.gv; g];
                obj.lbg = [obj.lbg; glb];
                obj.ubg = [obj.ubg; gub];
            end
        end
        function cost(obj)
            import casadi.*
            obj.Jk = obj.B(2:end)' * obj.Jj; % obj.B(r+1) dJ = obj.qr * obj.Zk(end)*sys.parameters.dt_scaling;
        end
        % Method to build map for constraint evaluation
        function build_map(obj, options)
            arguments
                obj; options.JN = obj.N;
            end
            import casadi.*
            obj.ng = length(obj.gv);
            obj.Eq_g_s = Function('Eq_g_s', {obj.x_1,obj.x, obj.xc,obj.u_1,obj.u, obj.uc, obj.z_1, obj.z, obj.zc, obj.p, obj.wg}, {obj.gv_s});
            obj.Eq_g_e = Function('Eq_g_e', {obj.x_1,obj.x, obj.xc,obj.u_1,obj.u, obj.uc, obj.z_1, obj.z, obj.zc, obj.p, obj.wg}, {obj.gv_e});
            obj.Eq_g_c = Function('Eq_g_c', {obj.x0, obj.x_end, obj.u0, obj.u_end, obj.z0, obj.z_end, obj.p, obj.wg}, {obj.gv_c});           
            obj.Eq_g = Function('Eq_g', {obj.x_1,obj.x, obj.xc, obj.u_1, obj.u, obj.uc, obj.z_1, obj.z, obj.zc, obj.pg, obj.p, obj.wg}, {obj.gv});
            obj.Eq_g_map = obj.Eq_g.map(obj.N, 'thread', obj.Nthread);
            obj.Eq_bg = Function('Eq_lbg', {obj.pgb}, {obj.lbg, obj.ubg});
            obj.Eq_bg_map = obj.Eq_bg.map(obj.N, 'thread', obj.Nthread); 
            if isempty(obj.Jk)
                obj.Jk = obj.Jj;
            end
            obj.J_fun = Function('J_fun', {obj.x_1,obj.x, obj.xc, obj.u_1, obj.u, obj.uc, obj.z_1, obj.z, obj.zc, obj.pg, obj.p, obj.wg}, {obj.Jk});
            %obj.J_fun_map = obj.J_fun.map(obj.N, 'thread', obj.Nthread);  
            obj.J_fun_map = obj.J_fun.map(options.JN, 'thread', obj.Nthread);  
            if strcmp(obj.sym_type, 'SXMX')
                obj.p = MX.sym('p', obj.np);
            end
        end
        function build_g(obj)
            import casadi.*
            obj.g = obj.Eq_g_map(obj.X_1,obj.X,obj.Xc,obj.U_1,obj.U, obj.Uc,obj.Z_1,obj.Z, obj.Zc, obj.Pg, obj.p, obj.Wg);
            obj.g = obj.g(:);
            [obj.lbg, obj.ubg] = obj.Eq_bg_map(obj.Pgb);
            obj.lbg = obj.lbg(:);
            obj.ubg = obj.ubg(:);
            if ~isempty(obj.gv_s)
                g_s = obj.Eq_g_s(obj.X_1(:,1),obj.X(:,1),obj.Xc(:,1:obj.d), obj.U_1(:,1), obj.U(:,1), obj.Uc(:,1:obj.d), obj.Z_1(:,1), obj.Z(:,1), obj.Zc(:,1:obj.d), obj.p, obj.Wg);
                obj.g = [g_s(:); obj.g];
                obj.lbg = [obj.lbg_s; obj.lbg];
                obj.ubg = [obj.ubg_s; obj.ubg];
            end
            if ~isempty(obj.gv_e)
                g_e = obj.Eq_g_e(obj.X_1(:,end),obj.X(:,end),obj.Xc(:,end-obj.d+1:end), obj.U_1(:,end), obj.U(:,end), obj.Uc(:,end-obj.d+1:end), obj.Z_1(:,end), obj.Z(:,end), obj.Zc(:,end-obj.d+1:end), obj.p, obj.Wg);
                obj.g = [obj.g; g_e(:)];
                obj.lbg = [obj.lbg; obj.lbg_e];
                obj.ubg = [obj.ubg; obj.ubg_e];
            end
            if ~isempty(obj.gv_c)
               g_c = obj.Eq_g_c(obj.X_1(:,1), obj.X(:,end), obj.U(:,1), obj.U(:,end), obj.Z(:,1), obj.Z(:,end), obj.p, obj.Wg); 
               obj.g = [obj.g; g_c(:)];
               obj.lbg = [obj.lbg; obj.lbg_c];
               obj.ubg = [obj.ubg; obj.ubg_c];
            end
            if ~isempty(obj.gv_m)
               g_m = obj.gv_m; 
               obj.g = [obj.g; g_m(:)];
               obj.lbg = [obj.lbg; obj.lbg_m];
               obj.ubg = [obj.ubg; obj.ubg_m];
            end            
        end
        function build_J(obj, options)
            arguments
                obj; options.range = 1:obj.N;
            end
            range_colloc = (options.range(1)*2 - 1):options.range(end)*2;
            import casadi.*
            %obj.J = sum(obj.J_fun_map(obj.X_1,obj.X,obj.Xc,obj.U_1,obj.U,obj.Uc,obj.Z_1,obj.Z,obj.Zc,obj.Pg,obj.p, obj.Wg));
            obj.J = sum(obj.J_fun_map(obj.X_1(:,options.range),obj.X(:,options.range),obj.Xc(:,range_colloc),obj.U_1(:,options.range),...
                obj.U(:,options.range),obj.Uc(:,range_colloc),obj.Z_1(:,options.range),obj.Z(:,options.range),...
                obj.Zc(:,range_colloc),obj.Pg(:,options.range),obj.p, obj.Wg));
        end
        % Method to build collocation coefficient
        function col_points(obj)
            % Collocation points definition
            import casadi.*
            if obj.d == 0
                obj.B = [1,1];
                obj.C = [];
                obj.D = [];
            else
                % Get collocation points
                tau_root = [0 collocation_points(obj.d, obj.colloc_type)];
                % Coefficients of the collocation equation
                obj.C = zeros(obj.d+1,obj.d+1);
                % Coefficients of the continuity equation
                obj.D = zeros(obj.d+1, 1);
                % Coefficients of the quadrature function
                obj.B = zeros(obj.d+1, 1);
                % Construct polynomial basis
                for j=1:obj.d+1
                    % Construct Lagrange polynomials to get the polynomial basis at the
                    % collocation point
                    coeff = 1;
                    for r=1:obj.d+1
                        if r ~= j
                            coeff = conv(coeff, [1, -tau_root(r)]);
                            coeff = coeff / (tau_root(j)-tau_root(r));
                        end
                    end
                    % Evaluate the polynomial at the final time to get the coefficients of
                    % the continuity equation
                    obj.D(j) = polyval(coeff, 1.0);
                    % Evaluate the time derivative of the polynomial at all collocation
                    % points to get the coefficients of the continuity equation
                    pder = polyder(coeff);
                    for r=1:obj.d+1
                        obj.C(j,r) = polyval(pder, tau_root(r));
                    end
                    % Evaluate the integral of the polynomial to get the coefficients of the
                    % quadrature function
                    pint = polyint(coeff);
                    obj.B(j) = polyval(pint, 1.0);
                end

                obj.xc_end = obj.D(1) * obj.x_1;
                for r = 1:obj.d
                    xpj{r} = obj.C(1,r+1) * obj.x_1;
                    for j = 1:obj.d
                        xpj{r} = xpj{r} + obj.C(j+1,r+1) * obj.xc(:,j);
                    end
                    obj.xc_end = obj.xc_end + obj.D(r+1) * obj.xc(:,r);
                end
                obj.xp = horzcat(xpj{:});
            end
        end
        function SS = empty_struct(obj, S, options)
            arguments
                obj; S; options.x0 = [];
            end
            if ~isfield(S,'lb')
                S.lb = [];
            end
            if ~isfield(S,'ub')
                S.ub = [];
            end            
            if isempty(options.x0) 
                SS = S;
            else
                if ~isfield(S, 'x0')
                    S.x0 = [];
                end                    
                SS = S;
            end
        end
        function casadi_sym(obj)
            import casadi.*
            if strcmp(obj.sym_type, 'SXMX')
                obj.w = MX.sym('w', (obj.nx+obj.nu+obj.nz+(obj.nx + obj.nuc + obj.nzc)*obj.d)*obj.N + obj.nx + obj.nwg, 1);
                sym_type_local = 'SX';
            else
                obj.w = casadi.(obj.sym_type).sym('w', (obj.nx+obj.nu+obj.nz+(obj.nx + obj.nuc + obj.nzc)*obj.d)*obj.N + obj.nx + obj.nwg, 1);
                sym_type_local = obj.sym_type;
            end
            obj.wg = casadi.(sym_type_local).sym('wg', obj.nwg); 
            obj.p = casadi.(sym_type_local).sym('p', obj.np);
            obj.x0 = casadi.(sym_type_local).sym('x0', obj.nx);
            obj.x_end = casadi.(sym_type_local).sym('x_end', obj.nx);
            obj.x = casadi.(sym_type_local).sym('x', obj.nx);
            obj.x_1 = casadi.(sym_type_local).sym('x_1', obj.nx);
            obj.xc = casadi.(sym_type_local).sym('xc', obj.nx, obj.d);
            obj.u = casadi.(sym_type_local).sym('u', obj.nu);
            obj.uc = casadi.(sym_type_local).sym('uc', obj.nuc, obj.d);
            obj.u0 = casadi.(sym_type_local).sym('u0', obj.nu);                    
            obj.u_end = casadi.(sym_type_local).sym('u_end', obj.nu);              
            obj.u_1 = casadi.(sym_type_local).sym('u_1', obj.nu);
            obj.z = casadi.(sym_type_local).sym('z', obj.nz);
            obj.zc = casadi.(sym_type_local).sym('zc', obj.nzc, obj.d);
            obj.z0 = casadi.(sym_type_local).sym('z0', obj.nz);           
            obj.z_end = casadi.(sym_type_local).sym('z_end', obj.nz);           
            obj.z_1 = casadi.(sym_type_local).sym('z_1', obj.nz);
            obj.pg = casadi.(sym_type_local).sym('pg', obj.Npg);
            obj.pgb = casadi.(sym_type_local).sym('pgb', obj.Npgb);
        end
        function build_solver(obj, opts, type)
            arguments
                obj; opts = []; type = [];
            end
            import casadi.*
            prob = struct('f', obj.J, 'x', obj.w, 'g', obj.g, 'p', obj.p);
            if isempty(opts)                
                if strcmp(type, 'integrator')
                    opts = struct;
                    opts.ipopt.linear_solver = 'ma57';
                    opts.ipopt.max_iter = 5000;
                    opts.ipopt.print_level = 5;
                    opts.ipopt.tol = 1e-4;
                elseif strcmp(type, 'NLP')
                    opts = struct( ...
                    'ipopt', ...
                    struct( ...
                    'max_iter', 5000, ... 
                    'linear_solver', 'ma57', ...
                    'print_level', 5, ...  
                    'nlp_scaling_method', 'none', ...
                    'ma57_pre_alloc', 2.5, ...
                    'mu_strategy', 'adaptive', ...
                    'tol', 1e-6), ...
                    'print_time', true);  
                elseif strcmp(type, 'default')
                    opts = struct();
                end
            end                
            obj.solver = nlpsol('solver', 'ipopt', prob, opts);            
        end
        function solution(obj, p)
            % Some other tests
            arguments
                obj; p = [];
            end
            if obj.np == 0
                obj.sol = obj.solver('x0', obj.w0, 'lbx', obj.lbw, 'ubx', obj.ubw, 'lbg', obj.lbg, 'ubg', obj.ubg);
            elseif strcmp(obj.check_x0, 'parametric')
                obj.sol = obj.solver('x0', obj.w0, 'lbx', obj.lbw(p), 'ubx', obj.ubw(p), 'lbg', obj.lbg, 'ubg', obj.ubg, 'p', p);
            else
                obj.sol = obj.solver('x0', obj.w0, 'lbx', obj.lbw, 'ubx', obj.ubw, 'lbg', obj.lbg, 'ubg', obj.ubg, 'p', p);
            end
            obj.X0_num = full(obj.sol.x(1:obj.nx));
            obj.sol_num = reshape(full(obj.sol.x(obj.nx+1:end-obj.nwg)), obj.nu + obj.nz + (obj.nuc + obj.nzc + obj.nx)*obj.d + obj.nx, obj.N);
            obj.w_global = full(obj.sol.x(end-obj.nwg + 1:end));
        end
    end
end

