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
            disp(obj.model_type)
        end        
        function aero_fun(obj)
            import casadi.*
            u = SX.sym('u');
            Fxaero = 0.5*obj.data.cx*obj.data.S*obj.data.rho*u^2;
            if isfield(obj.data, 'ab')
                cz1 = obj.data.ab*obj.data.cz;
                cz2 = (1-obj.data.ab)*obj.data.cz;
                Fzaero_1 = 0.5*cz1*obj.data.S*obj.data.rho*u^2;
                Fzaero_2 = 0.5*cz2*obj.data.S*obj.data.rho*u^2;
                obj.Faero = Function('Faero', {u}, {Fxaero, Fzaero_1, Fzaero_2}, {'u'}, {'Fx', 'Fz_1', 'Fz_2'});
            else
                Fzaero = 0.5*obj.data.cz*obj.data.S*obj.data.rho*u^2;
                obj.Faero = Function('Faero', {u}, {Fxaero, Fzaero}, {'u'}, {'Fx', 'Fz'});
            end
        end
        
    end
end