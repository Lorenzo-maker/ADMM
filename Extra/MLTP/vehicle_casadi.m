classdef vehicle_casadi < handle
    properties
        model_type
        data
        nx
        Xb
        nu
        Ub
        nuc
        nz
        Zb
        nzc
        nwg
        neq
        F 
        Faero
        Eq
        Q
        ndiseq
        fun_g3D
        L
    end
    
    methods
        function obj = vehicle_casadi(model_type, Data)
            arguments
                model_type = 'None';
                Data = [];
            end
            obj.model_type = model_type;
            obj.data = Data;
            switch model_type
                case 'point-mass'
                    obj.point_mass()
                case 'single-track'
                    obj.single_track()
                case 'double-track'
                    obj.double_track()
                case 'double-track-full'
                    obj.double_track_full()                    
                case 'double-track-spatial'
                    obj.double_track_spatial()
                otherwise
                    disp('no model selected')
            end
        end
        
        function point_mass(obj)
            import casadi.*         
            % Define global variables
            obj.nwg = 0;
            % Define vector state [u, v, r, x, y, psi]
            obj.nx = 6;
            X = SX.sym('X', obj.nx, 1);
            X_un = X.*obj.data.X_scale;
            obj.Xb.lb = [0; -1; -1; -1; -1; -1];            
            obj.Xb.ub = [1; 1; 1; 1; 1; 1];
            % Define controls vector [Fx, rp]
            obj.nu = 2;
            U = SX.sym('U', obj.nu, 1);
            U_un = U.*obj.data.U_scale;
            obj.Ub.lb = [-1; -1];          
            obj.Ub.ub = [1; 1];
            % Define controls in collocations
            obj.nuc = 0;
            % Define algebraic parameters [Fz, Fy, Pow, ep, h]
            obj.nz = 5;
            Z = SX.sym('Z', obj.nz, 1);
            Z_un = Z.*obj.data.Z_scale;
            obj.Zb.lb = [0; -1; obj.data.Pmin/obj.data.P_scale; -1; 0];           
            obj.Zb.ub = [1; 1; obj.data.Pmax/obj.data.P_scale; 1; 1];
            % Define algebraic in collocations
            obj.nzc = 0;
            % Aero function
            obj.aero_fun();
            F_aero = obj.Faero('u', X_un(1));
            % Define ODE
            u_dot = (U_un(1) - F_aero.Fx)/obj.data.m + X_un(2)*X_un(3);
            v_dot = Z_un(2)/obj.data.m - X_un(1)*X_un(3);
            r_dot = U_un(2);
            x_dot = ((X_un(1))*cos(X_un(6)) - (X_un(2))*sin(X_un(6)));
            y_dot = ((X_un(1))*sin(X_un(6)) + (X_un(2))*cos(X_un(6)));
            psi_dot = X_un(3);
            X_dot = [u_dot;v_dot;r_dot;x_dot;y_dot;psi_dot];
            obj.F = Function('F', {X, U, Z}, {X_dot});
            % Define Eq algebraic
            eq1 = (Z_un(1) - (obj.data.m*obj.data.G + F_aero.Fz))/obj.data.Z_scale(1);
            eq2 = (U_un(1)*X_un(1) - Z_un(3))./obj.data.Z_scale(3);
            eq  = [eq1; eq2];
            obj.Eq = Function('Eq', {X, U, Z}, {eq});
            obj.neq   = length(eq);            
        end
        function single_track(obj)
            disp(obj.model_type)
        end
        function double_track(obj)
            import casadi.*         
            % Define global variables
            obj.nwg = 0;
            % Define vector state [u, v, r, x, y, psi, omega_ij]
            obj.nx = 6 + 4;
            X = SX.sym('X', obj.nx, 1);
            X_un = X.*obj.data.X_scale;
            obj.Xb.lb = [0; -1; -1; -1; -1; -1; 0; 0; 0; 0];            
            obj.Xb.ub = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1];
            % Define controls vector [Ta, Tb, delta]
            obj.nu = 3;
            U = SX.sym('U', obj.nu, 1);
            U_un = U.*obj.data.U_scale;
            obj.Ub.lb = [0; -1; -1];          
            obj.Ub.ub = [1; 0; 1];            
            % Define controls in collocations
            obj.nuc = 0;
            % Define algebraic parameters [Fz11, Fz12, Fz21, Fz22, Fx11, Fx12, Fx21, Fx22, Fy11, Fy12, Fy21, Fy22, ep, h]
            obj.nz = 6;
            Z = SX.sym('Z', obj.nz, 1);
            Z_un = Z.*obj.data.Z_scale;
            obj.Zb.lb = [0; 0; 0; 0; -1; 1e-2];           
            obj.Zb.ub = [1; 1; 1; 1; 1; 1];
            % Define algebraic in collocations
            obj.nzc = 0;
            % Aero function
            obj.aero_fun();
            F_aero = obj.Faero('u', X_un(1));
            % Define 3D gravity
            g3D = SX.sym('g3D', 3, 1);
            Rotsg = SX.sym('Rotsg',3,3);
            g3D_sym   = [cos(X_un(6)), sin(X_un(6)), 0; -sin(X_un(6)), cos(X_un(6)), 0; 0, 0, 1]*Rotsg*[0,0,obj.data.G]';
            obj.fun_g3D = Function('g3D',{X(6), Rotsg}, {g3D_sym}); 
            % Define Steering law
            d11 = -obj.data.d0 + obj.data.tau*U_un(3) + obj.data.ep1*obj.data.t1*(obj.data.tau*U_un(3))^2/(2*obj.data.l);
            d12 = obj.data.d0 + obj.data.tau*U_un(3) - obj.data.ep1*obj.data.t1*(obj.data.tau*U_un(3))^2/(2*obj.data.l);
            % Define Vxij and Vyij
            Fx_tyre = obj.data.tyre.Fx;
            Fy_tyre = obj.data.tyre.Fy;
            Vx11 = (X_un(1) - X_un(3)*obj.data.t1/2)*cos(d11) + (X_un(2) + X_un(3)*obj.data.a1)*sin(d11);
            Vx12 = (X_un(1) + X_un(3)*obj.data.t1/2)*cos(d12) + (X_un(2) + X_un(3)*obj.data.a1)*sin(d12);
            Vx21 = X_un(1) - X_un(3)*obj.data.t2/2;
            Vx22 = X_un(1) + X_un(3)*obj.data.t2/2;
            Vy11 = (X_un(2) + X_un(3)*obj.data.a1)*cos(d11) - (X_un(1) - X_un(3)*obj.data.t1/2)*sin(d11);
            Vy12 = (X_un(2) + X_un(3)*obj.data.a1)*cos(d12) - (X_un(1) + X_un(3)*obj.data.t1/2)*sin(d12);
            Vy21 = X_un(2) - X_un(3)*obj.data.a2;
            Vy22 = X_un(2) - X_un(3)*obj.data.a2;
            %Slip k
            k11 =  1 - X_un(7)*obj.data.rw/(Vx11 + 1e-3); % 1 - (X_un(7)*obj.data.rw)/(X_un(1)+1e-3);
            k12 =  1 - X_un(8)*obj.data.rw/(Vx12 + 1e-3); %1 - (X_un(8)*obj.data.rw)/(X_un(1)+1e-3);
            k21 =  1 - X_un(9)*obj.data.rw/(Vx21 + 1e-3); %1 - (X_un(9)*obj.data.rw)/(X_un(1)+1e-3);
            k22 =  1 - X_un(10)*obj.data.rw/(Vx22 + 1e-3); %1 - (X_un(10)*obj.data.rw)/(X_un(1)+1e-3);
            %slip alpha
            alpha11 = atan(-Vy11/(Vx11 + 1e-3)); %atan(-((X_un(2) + obj.data.a1*X_un(3))/(X_un(1)+1e-3) + tan(U_un(3)))); 
            alpha12 = atan(-Vy12/(Vx12 + 1e-3)); %atan(-((X_un(2) + obj.data.a1*X_un(3))/(X_un(1)+1e-3) + tan(U_un(3)))); 
            alpha21 = atan(-Vy21/(Vx21 + 1e-3)); %atan(-((X_un(2) - obj.data.a2*X_un(3))/(X_un(1)+1e-3))); 
            alpha22 = atan(-Vy22/(Vx22 + 1e-3)); %atan(-((X_un(2) - obj.data.a2*X_un(3))/(X_un(1)+1e-3))); 
            % Define ODE  
            k = SX.sym('k');
            alph = SX.sym('alph');
            Fz = SX.sym('Fz');
            T = SX.sym('T');
            Fx = if_else_smooth(abs(T),...                                                % expression to check
                                0,...                                                % greater than
                                min([abs(Fx_tyre(k, alph, Fz, 0)), (abs(T)/obj.data.rw), 1.2.*Fz]).*sign(-k),... % true
                                Fx_tyre(k, alph, Fz, 0),...                                               % false
                                'C', 20);                                                     % smooth parameter
            Fx_fun = Function('Fx', {k, alph, Fz, T}, {Fx});  
            Fx11 = Fx_tyre(k11, alpha11, Z_un(1), 0);
            Fx12 = Fx_tyre(k12, alpha12, Z_un(2), 0);
            Fx21 = Fx_tyre(k21, alpha21, Z_un(3), 0);
            Fx22 = Fx_tyre(k22, alpha22, Z_un(4), 0);
            Fy11 = Fy_tyre(alpha11, k11, Z_un(1), 0);           
            Fy12 = Fy_tyre(alpha12, k12, Z_un(2), 0);           
            Fy21 = Fy_tyre(alpha21, k21, Z_un(3), 0);
            Fy22 = Fy_tyre(alpha22, k22, Z_un(4), 0);
            Fx_tot = Fx21 + Fx22 + (Fx11*cos(d11) + Fx12*cos(d12)) - (Fy11*sin(d11) + Fy12*sin(d12));
            Fy_tot = Fy21 + Fy22 + (Fy11*cos(d11) + Fy12*cos(d12)) + (Fx11*sin(d11) + Fx12*sin(d12));
            Mz_tot = (Fx22 - Fx21)*obj.data.t2/2 + (Fx12*cos(d12) - Fx11*cos(d11))*obj.data.t1/2 -(Fy21 + Fy22)*obj.data.a2 + ((Fy11*cos(d11) + Fy12*cos(d12)) + (Fx11*sin(d11) + Fx12*sin(d12)))*obj.data.a1;                       
            % ODE
            etas = obj.data.eta^(if_else_smooth(X_un(10) - X_un(9), 0, 1, -1, 'C', 20));
            Tr = U_un(1)*etas/(1 + etas);
            Tl = U_un(1)/(1 + etas);
            u_dot = (Fx_tot - F_aero.Fx)/obj.data.m + X_un(2)*X_un(3) - g3D(1);
            v_dot = Fy_tot/obj.data.m - X_un(1)*X_un(3) - g3D(2);
            r_dot = Mz_tot/obj.data.Jz;
            x_dot = ((X_un(1))*cos(X_un(6)) - (X_un(2))*sin(X_un(6)));
            y_dot = ((X_un(1))*sin(X_un(6)) + (X_un(2))*cos(X_un(6)));
            psi_dot = X_un(3);
            omega11_dot = (obj.data.kb*U_un(2)/2 - Fx11*obj.data.rw)/obj.data.Jw;
            omega12_dot = (obj.data.kb*U_un(2)/2 - Fx12*obj.data.rw)/obj.data.Jw;
            omega21_dot = (Tl + (1-obj.data.kb)*U_un(2)/2 - Fx21*obj.data.rw)/obj.data.Jw;
            omega22_dot = (Tr + (1-obj.data.kb)*U_un(2)/2 - Fx22*obj.data.rw)/obj.data.Jw;
            X_dot = [u_dot;v_dot;r_dot;x_dot;y_dot;psi_dot;omega11_dot;omega12_dot;omega21_dot;omega22_dot];
            obj.F = Function('F', {X, U, Z, g3D}, {X_dot});           
            % Eq
            Dz = obj.data.m*(u_dot - X_un(2)*X_un(3))*obj.data.hg/obj.data.l;
            Y = obj.data.m*(v_dot + X_un(1)*X_un(3));
            Y1 = (Fy11*cos(d11) + Fy12*cos(d12)) + (Fx11*sin(d11) + Fx12*sin(d12));
            Y2 = Fy21 + Fy22;
            Dz1 = (1/obj.data.t1)*(obj.data.k_phi1*(Y*(obj.data.hg - obj.data.qm))/obj.data.k_phi + Y1*obj.data.q1);
            Dz2 = (1/obj.data.t2)*(obj.data.k_phi2*(Y*(obj.data.hg - obj.data.qm))/obj.data.k_phi + Y2*obj.data.q2);           
            eq1 = (Z_un(1) - (0.5*(obj.data.m*g3D(3)*obj.data.a2/obj.data.l + F_aero.Fz_1 - Dz) - Dz1))/obj.data.Z_scale(1);
            eq2 = (Z_un(2) - (0.5*(obj.data.m*g3D(3)*obj.data.a2/obj.data.l + F_aero.Fz_1 - Dz) + Dz1))/obj.data.Z_scale(2);
            eq3 = (Z_un(3) - (0.5*(obj.data.m*g3D(3)*obj.data.a1/obj.data.l + F_aero.Fz_2 + Dz) - Dz2))/obj.data.Z_scale(3);
            eq4 = (Z_un(4) - (0.5*(obj.data.m*g3D(3)*obj.data.a1/obj.data.l + F_aero.Fz_2 + Dz) + Dz2))/obj.data.Z_scale(4);   
            eq  = [eq1; eq2; eq3; eq4];
            obj.Eq = Function('Eq', {X, U, Z, g3D}, {eq});
            obj.neq   = length(eq);
            % quantities
            q1 = k11; q2 = k12; q3 = k21; q4 = k22;
            q5 = alpha11; q6 = alpha12; q7 = alpha21; q8 = alpha22;
            q = [q1; q2; q3; q4; q5; q6; q7; q8];
            obj.Q = Function('Q', {X, U, Z}, {q});
        end     
        function double_track_full(obj)
            import casadi.*         
            % Define global variables
            obj.nwg = 0;
            % Define vector state [u, v, r, x, y, psi, omega_ij]
            obj.nx = 6 + 4;
            X = SX.sym('X', obj.nx, 1);
            X_un = X.*obj.data.X_scale;
            obj.Xb.lb = [0; -1; -1; -1; -1; -1; 0; 0; 0; 0];            
            obj.Xb.ub = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1];
            % Define controls vector [Ta, Tb, delta]
            obj.nu = 3;
            U = SX.sym('U', obj.nu, 1);
            U_un = U.*obj.data.U_scale;
            obj.Ub.lb = [0; -1; -1];          
            obj.Ub.ub = [1; 0; 1];            
            % Define controls in collocations
            obj.nuc = 0;
            % Define algebraic parameters [Fz11, Fz12, Fz21, Fz22, Fx11, Fx12, Fx21, Fx22, Fy11, Fy12, Fy21, Fy22, ep, h]
            obj.nz = 14;
            Z = SX.sym('Z', obj.nz, 1);
            Z_un = Z.*obj.data.Z_scale;
            obj.Zb.lb = [0;0;0;0;-1*ones(8,1);-1; 1e-2];           
            obj.Zb.ub = [1*ones(12,1);1; 1];
            % Define algebraic in collocations
            obj.nzc = 0;
            % Aero function
            obj.aero_fun();
            F_aero = obj.Faero('u', X_un(1));
            % Define 3D gravity
            g3D = SX.sym('g3D', 3, 1);
            Rotsg = SX.sym('Rotsg',3,3);
            g3D_sym   = [cos(X_un(6)), sin(X_un(6)), 0; -sin(X_un(6)), cos(X_un(6)), 0; 0, 0, 1]*Rotsg*[0,0,obj.data.G]';
            obj.fun_g3D = Function('g3D',{X(6), Rotsg}, {g3D_sym}); 
            % Define Steering law
            d11 = -obj.data.d0 + obj.data.tau*U_un(3) + obj.data.ep1*obj.data.t1*(obj.data.tau*U_un(3))^2/(2*obj.data.l);
            d12 = obj.data.d0 + obj.data.tau*U_un(3) - obj.data.ep1*obj.data.t1*(obj.data.tau*U_un(3))^2/(2*obj.data.l);
            % Define Vxij and Vyij
            Fx_tyre = obj.data.tyre.Fx;
            Fy_tyre = obj.data.tyre.Fy;
            Vx11 = (X_un(1) - X_un(3)*obj.data.t1/2)*cos(d11) + (X_un(2) + X_un(3)*obj.data.a1)*sin(d11);
            Vx12 = (X_un(1) + X_un(3)*obj.data.t1/2)*cos(d12) + (X_un(2) + X_un(3)*obj.data.a1)*sin(d12);
            Vx21 = X_un(1) - X_un(3)*obj.data.t2/2;
            Vx22 = X_un(1) + X_un(3)*obj.data.t2/2;
            Vy11 = (X_un(2) + X_un(3)*obj.data.a1)*cos(d11) - (X_un(1) - X_un(3)*obj.data.t1/2)*sin(d11);
            Vy12 = (X_un(2) + X_un(3)*obj.data.a1)*cos(d12) - (X_un(1) + X_un(3)*obj.data.t1/2)*sin(d12);
            Vy21 = X_un(2) - X_un(3)*obj.data.a2;
            Vy22 = X_un(2) - X_un(3)*obj.data.a2;
            %Slip k
            k11 =  1 - X_un(7)*obj.data.rw/(Vx11 + 1e-3); % 1 - (X_un(7)*obj.data.rw)/(X_un(1)+1e-3);
            k12 =  1 - X_un(8)*obj.data.rw/(Vx12 + 1e-3); %1 - (X_un(8)*obj.data.rw)/(X_un(1)+1e-3);
            k21 =  1 - X_un(9)*obj.data.rw/(Vx21 + 1e-3); %1 - (X_un(9)*obj.data.rw)/(X_un(1)+1e-3);
            k22 =  1 - X_un(10)*obj.data.rw/(Vx22 + 1e-3); %1 - (X_un(10)*obj.data.rw)/(X_un(1)+1e-3);
            %slip alpha
            alpha11 = atan(-Vy11/(Vx11 + 1e-3)); %atan(-((X_un(2) + obj.data.a1*X_un(3))/(X_un(1)+1e-3) + tan(U_un(3)))); 
            alpha12 = atan(-Vy12/(Vx12 + 1e-3)); %atan(-((X_un(2) + obj.data.a1*X_un(3))/(X_un(1)+1e-3) + tan(U_un(3)))); 
            alpha21 = atan(-Vy21/(Vx21 + 1e-3)); %atan(-((X_un(2) - obj.data.a2*X_un(3))/(X_un(1)+1e-3))); 
            alpha22 = atan(-Vy22/(Vx22 + 1e-3)); %atan(-((X_un(2) - obj.data.a2*X_un(3))/(X_un(1)+1e-3))); 
            % Define ODE  
            Fx11 = Z_un(5);%Fx_tyre(k11, alpha11, Z_un(1), 0);
            Fx12 = Z_un(6);%Fx_tyre(k12, alpha12, Z_un(2), 0);
            Fx21 = Z_un(7);%Fx_tyre(k21, alpha21, Z_un(3), 0);
            Fx22 = Z_un(8);%Fx_tyre(k22, alpha22, Z_un(4), 0);
            Fy11 = Z_un(9);%Fy_tyre(alpha11, k11, Z_un(1), 0);           
            Fy12 = Z_un(10);%Fy_tyre(alpha12, k12, Z_un(2), 0);           
            Fy21 = Z_un(11);%Fy_tyre(alpha21, k21, Z_un(3), 0);
            Fy22 = Z_un(12);%Fy_tyre(alpha22, k22, Z_un(4), 0);
            Fx_tot = Fx21 + Fx22 + (Fx11*cos(d11) + Fx12*cos(d12)) - (Fy11*sin(d11) + Fy12*sin(d12));
            Fy_tot = Fy21 + Fy22 + (Fy11*cos(d11) + Fy12*cos(d12)) + (Fx11*sin(d11) + Fx12*sin(d12));
            Mz_tot = (Fx22 - Fx21)*obj.data.t2/2 + (Fx12*cos(d12) - Fx11*cos(d11))*obj.data.t1/2 -(Fy21 + Fy22)*obj.data.a2 + ((Fy11*cos(d11) + Fy12*cos(d12)) + (Fx11*sin(d11) + Fx12*sin(d12)))*obj.data.a1;                       
            % ODE
            etas = obj.data.eta^(if_else_smooth(X_un(10) - X_un(9), 0, 1, -1, 'C', 20));
            Tr = U_un(1)*etas/(1 + etas);
            Tl = U_un(1)/(1 + etas);
            u_dot = (Fx_tot - F_aero.Fx)/obj.data.m + X_un(2)*X_un(3) - g3D(1);
            v_dot = Fy_tot/obj.data.m - X_un(1)*X_un(3) - g3D(2);
            r_dot = Mz_tot/obj.data.Jz;
            x_dot = ((X_un(1))*cos(X_un(6)) - (X_un(2))*sin(X_un(6)));
            y_dot = ((X_un(1))*sin(X_un(6)) + (X_un(2))*cos(X_un(6)));
            psi_dot = X_un(3);
            omega11_dot = (obj.data.kb*U_un(2)/2 - Fx11*obj.data.rw)/obj.data.Jw;
            omega12_dot = (obj.data.kb*U_un(2)/2 - Fx12*obj.data.rw)/obj.data.Jw;
            omega21_dot = (Tl + (1-obj.data.kb)*U_un(2)/2 - Fx21*obj.data.rw)/obj.data.Jw;
            omega22_dot = (Tr + (1-obj.data.kb)*U_un(2)/2 - Fx22*obj.data.rw)/obj.data.Jw;
            X_dot = [u_dot;v_dot;r_dot;x_dot;y_dot;psi_dot;omega11_dot;omega12_dot;omega21_dot;omega22_dot];
            obj.F = Function('F', {X, U, Z, g3D}, {X_dot});           
            % Eq
            Dz = obj.data.m*(u_dot - X_un(2)*X_un(3))*obj.data.hg/obj.data.l;
            Y = obj.data.m*(v_dot + X_un(1)*X_un(3));
            Y1 = (Fy11*cos(d11) + Fy12*cos(d12)) + (Fx11*sin(d11) + Fx12*sin(d12));
            Y2 = Fy21 + Fy22;
            Dz1 = (1/obj.data.t1)*(obj.data.k_phi1*(Y*(obj.data.hg - obj.data.qm))/obj.data.k_phi + Y1*obj.data.q1);
            Dz2 = (1/obj.data.t2)*(obj.data.k_phi2*(Y*(obj.data.hg - obj.data.qm))/obj.data.k_phi + Y2*obj.data.q2);           
            eq1 = (Z_un(1) - (0.5*(obj.data.m*g3D(3)*obj.data.a2/obj.data.l + F_aero.Fz_1 - Dz) - Dz1))/obj.data.Z_scale(1);
            eq2 = (Z_un(2) - (0.5*(obj.data.m*g3D(3)*obj.data.a2/obj.data.l + F_aero.Fz_1 - Dz) + Dz1))/obj.data.Z_scale(2);
            eq3 = (Z_un(3) - (0.5*(obj.data.m*g3D(3)*obj.data.a1/obj.data.l + F_aero.Fz_2 + Dz) - Dz2))/obj.data.Z_scale(3);
            eq4 = (Z_un(4) - (0.5*(obj.data.m*g3D(3)*obj.data.a1/obj.data.l + F_aero.Fz_2 + Dz) + Dz2))/obj.data.Z_scale(4); 
            eq5 =  Z_un(5) - Fx_tyre(k11, alpha11, Z_un(1), 0);
            eq6 =  Z_un(6) - Fx_tyre(k12, alpha12, Z_un(2), 0);
            eq7 =  Z_un(7) - Fx_tyre(k21, alpha21, Z_un(3), 0);
            eq8 =  Z_un(8) - Fx_tyre(k22, alpha22, Z_un(4), 0);
            eq9 =  Z_un(9) - Fy_tyre(alpha11, k11, Z_un(1), 0);           
            eq10 = Z_un(10)- Fy_tyre(alpha12, k12, Z_un(2), 0);           
            eq11 = Z_un(11)- Fy_tyre(alpha21, k21, Z_un(3), 0);
            eq12 = Z_un(12)- Fy_tyre(alpha22, k22, Z_un(4), 0); 
            eq  = [eq1; eq2; eq3; eq4; eq5; eq6; eq7; eq8; eq9; eq10; eq11; eq12];
            obj.Eq = Function('Eq', {X, U, Z, g3D}, {eq});
            obj.neq   = length(eq);
            % quantities
            q1 = k11; q2 = k12; q3 = k21; q4 = k22;
            q5 = alpha11; q6 = alpha12; q7 = alpha21; q8 = alpha22;
            q = [q1; q2; q3; q4; q5; q6; q7; q8];
            obj.Q = Function('Q', {X, U, Z}, {q});
        end             
        function double_track_spatial(obj)
            import casadi.*         
            % Define global variables
            obj.nwg = 0;
            % Define vector state [u, v, r, ep, epsi, omega_ij]
            obj.nx = 5 + 4;
            X = SX.sym('X', obj.nx, 1);
            X_un = X.*obj.data.X_scale;
            obj.Xb.lb = [0; -1; -1; -1; -1; 0; 0; 0; 0];            
            obj.Xb.ub = [1; 1; 1; 1; 1; 1; 1; 1; 1];
            % Define controls vector [Ta, Tb, delta]
            obj.nu = 3;
            U = SX.sym('U', obj.nu, 1);
            U_un = U.*obj.data.U_scale;
            obj.Ub.lb = [0; -1; -1];          
            obj.Ub.ub = [1; 0; 1];
            % Define controls in collocations
            obj.nuc = 0;
            % Define algebraic parameters [Fz11, Fz12, Fz21, Fz22]
            obj.nz = 4;
            Z = SX.sym('Z', obj.nz, 1);
            Z_un = Z.*obj.data.Z_scale;
            obj.Zb.lb = [0; 0; 0; 0];           
            obj.Zb.ub = [1; 1; 1; 1];
            % Define algebraic in collocations
            obj.nzc = 0;
            % Aero function
            obj.aero_fun();
            F_aero = obj.Faero('u', X_un(1));
            % Define 3D gravity
            tol = 1e-4;
            g3D = SX.sym('g3D', 3, 1);
            Rotsg = SX.sym('Rotsg',3,3);
            g3D_sym   = [cos(X_un(5)+tol), sin(X_un(5)+tol), 0; -sin(X_un(5)+tol), cos(X_un(5)+tol), 0; 0, 0, 1]*Rotsg*[0,0,obj.data.G]';
            obj.fun_g3D = Function('g3D',{X(5), Rotsg}, {g3D_sym});  
            % Define steering law
            d11 = -obj.data.d0 + obj.data.tau*U_un(3) + obj.data.ep1*obj.data.t1*(obj.data.tau*U_un(3))^2/(2*obj.data.l);
            d12 = obj.data.d0 + obj.data.tau*U_un(3) - obj.data.ep1*obj.data.t1*(obj.data.tau*U_un(3))^2/(2*obj.data.l);
            % Define Vxij and Vyij
            Fx_tyre = obj.data.tyre.Fx;
            Fy_tyre = obj.data.tyre.Fy;
            Vx11 = (X_un(1) - X_un(3)*obj.data.t1/2)*cos(d11) + (X_un(2) + X_un(3)*obj.data.a1)*sin(d11);
            Vx12 = (X_un(1) + X_un(3)*obj.data.t1/2)*cos(d12) + (X_un(2) + X_un(3)*obj.data.a1)*sin(d12);
            Vx21 = X_un(1) - X_un(3)*obj.data.t2/2;
            Vx22 = X_un(1) + X_un(3)*obj.data.t2/2;
            Vy11 = (X_un(2) + X_un(3)*obj.data.a1)*cos(d11) - (X_un(1) - X_un(3)*obj.data.t1/2)*sin(d11);
            Vy12 = (X_un(2) + X_un(3)*obj.data.a1)*cos(d12) - (X_un(1) + X_un(3)*obj.data.t1/2)*sin(d12);
            Vy21 = X_un(2) - X_un(3)*obj.data.a2;
            Vy22 = X_un(2) - X_un(3)*obj.data.a2;
            %Slip k
            k11 =  1 - X_un(6)*obj.data.rw/(Vx11 + 1e-3); % 1 - (X_un(7)*obj.data.rw)/(X_un(1)+1e-3);
            k12 =  1 - X_un(7)*obj.data.rw/(Vx12 + 1e-3); %1 - (X_un(8)*obj.data.rw)/(X_un(1)+1e-3);
            k21 =  1 - X_un(8)*obj.data.rw/(Vx21 + 1e-3); %1 - (X_un(9)*obj.data.rw)/(X_un(1)+1e-3);
            k22 =  1 - X_un(9)*obj.data.rw/(Vx22 + 1e-3); %1 - (X_un(10)*obj.data.rw)/(X_un(1)+1e-3);
            %slip alpha
            alpha11 = atan(-Vy11/(Vx11 + 1e-3)); %atan(-((X_un(2) + obj.data.a1*X_un(3))/(X_un(1)+1e-3) + tan(U_un(3)))); 
            alpha12 = atan(-Vy12/(Vx12 + 1e-3)); %atan(-((X_un(2) + obj.data.a1*X_un(3))/(X_un(1)+1e-3) + tan(U_un(3)))); 
            alpha21 = atan(-Vy21/(Vx21 + 1e-3)); %atan(-((X_un(2) - obj.data.a2*X_un(3))/(X_un(1)+1e-3))); 
            alpha22 = atan(-Vy22/(Vx22 + 1e-3)); %atan(-((X_un(2) - obj.data.a2*X_un(3))/(X_un(1)+1e-3))); 
            % Define ODE  
            k = SX.sym('k');
            alph = SX.sym('alph');
            Fz = SX.sym('Fz');
            T = SX.sym('T');
            Fx = if_else_smooth(abs(T),...                                                % expression to check
                                0,...                                                % greater than
                                min([abs(Fx_tyre(k, alph, Fz, 0)), (abs(T)/obj.data.rw), 1.2.*Fz]).*sign(-k),... % true
                                Fx_tyre(k, alph, Fz, 0),...                                               % false
                                'C', 20);                                                     % smooth parameter
            Fx_fun = Function('Fx', {k, alph, Fz, T}, {Fx});  
            Fx11 = Fx_tyre(k11, alpha11, Z_un(1), 0);
            Fx12 = Fx_tyre(k12, alpha12, Z_un(2), 0);
            Fx21 = Fx_tyre(k21, alpha21, Z_un(3), 0);
            Fx22 = Fx_tyre(k22, alpha22, Z_un(4), 0);
            Fy11 = Fy_tyre(alpha11, k11, Z_un(1), 0);           
            Fy12 = Fy_tyre(alpha12, k12, Z_un(2), 0);           
            Fy21 = Fy_tyre(alpha21, k21, Z_un(3), 0);
            Fy22 = Fy_tyre(alpha22, k22, Z_un(4), 0);
            Fx_tot = Fx21 + Fx22 + (Fx11*cos(d11) + Fx12*cos(d12)) - (Fy11*sin(d11) + Fy12*sin(d12));
            Fy_tot = Fy21 + Fy22 + (Fy11*cos(d11) + Fy12*cos(d12)) + (Fx11*sin(d11) + Fx12*sin(d12));
            Mz_tot = (Fx22 - Fx21)*obj.data.t2/2 + (Fx12*cos(d12) - Fx11*cos(d11))*obj.data.t1/2 -(Fy21 + Fy22)*obj.data.a2 + ((Fy11*cos(d11) + Fy12*cos(d12)) + (Fx11*sin(d11) + Fx12*sin(d12)))*obj.data.a1;                       
            % sp
            kn = SX.sym('kn');
            normts = SX.sym('normts');
            sp    = (X_un(1)*cos(X_un(5))-(X_un(2))*sin(X_un(5)))/(1 - kn*(X_un(4)));
            fdtda = normts/(sp + tol);
            % ODE
            etas = obj.data.eta^(if_else_smooth(X_un(9) - X_un(8), 0, 1, -1, 'C', 20));
            Tr = U_un(1)*etas/(1 + etas);
            Tl = U_un(1)/(1 + etas);
            u_dot = ((Fx_tot - F_aero.Fx)/obj.data.m + X_un(2)*X_un(3) - g3D(1))*fdtda;
            v_dot = (Fy_tot/obj.data.m - X_un(1)*X_un(3) - g3D(2))*fdtda;
            r_dot = (Mz_tot/obj.data.Jz)*fdtda;           
            ep_dot = ((X_un(1))*sin(X_un(5)) + (X_un(2))*cos(X_un(5)))*fdtda;
            ef_dot = (X_un(3) - kn*sp)*fdtda;
            omega11_dot = ((obj.data.kb*U_un(2)/2 - Fx11*obj.data.rw)/obj.data.Jw)*fdtda;
            omega12_dot = ((obj.data.kb*U_un(2)/2 - Fx12*obj.data.rw)/obj.data.Jw)*fdtda;
            omega21_dot = ((Tl + (1-obj.data.kb)*U_un(2)/2 - Fx21*obj.data.rw)/obj.data.Jw)*fdtda;
            omega22_dot = ((Tr + (1-obj.data.kb)*U_un(2)/2 - Fx22*obj.data.rw)/obj.data.Jw)*fdtda;
            X_dot = [u_dot;v_dot;r_dot;ep_dot;ef_dot;omega11_dot;omega12_dot;omega21_dot;omega22_dot];
            obj.F = Function('F', {X, U, Z, g3D, kn, normts}, {X_dot});  
            % Cost function
            obj.L = Function('L', {X, kn, normts}, {fdtda});
            % Eq
            Dz = obj.data.m*(u_dot/fdtda - X_un(2)*X_un(3))*obj.data.hg/obj.data.l;
            Y = obj.data.m*(v_dot/fdtda + X_un(1)*X_un(3));
            Y1 = (Fy11*cos(d11) + Fy12*cos(d12)) + (Fx11*sin(d11) + Fx12*sin(d12));
            Y2 = Fy21 + Fy22;
            Dz1 = (1/obj.data.t1)*(obj.data.k_phi1*(Y*(obj.data.hg - obj.data.qm))/obj.data.k_phi + Y1*obj.data.q1);
            Dz2 = (1/obj.data.t2)*(obj.data.k_phi2*(Y*(obj.data.hg - obj.data.qm))/obj.data.k_phi + Y2*obj.data.q2);           
            eq1 = (Z_un(1) - (0.5*(obj.data.m*g3D(3)*obj.data.a2/obj.data.l + F_aero.Fz_1 - Dz) - Dz1))/obj.data.Z_scale(1);
            eq2 = (Z_un(2) - (0.5*(obj.data.m*g3D(3)*obj.data.a2/obj.data.l + F_aero.Fz_1 - Dz) + Dz1))/obj.data.Z_scale(2);
            eq3 = (Z_un(3) - (0.5*(obj.data.m*g3D(3)*obj.data.a1/obj.data.l + F_aero.Fz_2 + Dz) - Dz2))/obj.data.Z_scale(3);
            eq4 = (Z_un(4) - (0.5*(obj.data.m*g3D(3)*obj.data.a1/obj.data.l + F_aero.Fz_2 + Dz) + Dz2))/obj.data.Z_scale(4);   
            eq  = [eq1; eq2; eq3; eq4];
            obj.Eq = Function('Eq', {X, U, Z, g3D, kn, normts}, {eq});
            obj.neq   = length(eq);       
            % quantities
            q1 = k11; q2 = k12; q3 = k21; q4 = k22;
            q5 = alpha11; q6 = alpha12; q7 = alpha21; q8 = alpha22;
            q = [q1; q2; q3; q4; q5; q6; q7; q8];
            obj.Q = Function('Q', {X, U, Z}, {q});            
        end
        function aero_fun(obj)
            import casadi.*
            u = SX.sym('u');
            Fxaero = 0.5*obj.data.cx*obj.data.S*obj.data.rho*u^2;
            if strcmp(obj.model_type, 'point-mass')
                Fzaero = 0.5*obj.data.cz*obj.data.S*obj.data.rho*u^2;
                obj.Faero = Function('Faero', {u}, {Fxaero, Fzaero}, {'u'}, {'Fx', 'Fz'});                
            else 
                cz1 = obj.data.ab*obj.data.cz;
                cz2 = (1-obj.data.ab)*obj.data.cz;
                Fzaero_1 = 0.5*cz1*obj.data.S*obj.data.rho*u^2;
                Fzaero_2 = 0.5*cz2*obj.data.S*obj.data.rho*u^2;
                obj.Faero = Function('Faero', {u}, {Fxaero, Fzaero_1, Fzaero_2}, {'u'}, {'Fx', 'Fz_1', 'Fz_2'});
            end
        end

        
    end
end