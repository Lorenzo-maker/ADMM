classdef FSAE < handle

    
    properties         
        model_type
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
        sym
        data      
        dyn_fun
        out_fun
    end
    
    methods        
        function self = FSAE(track, alfarange, model_type)            
            
            % Set the default symbolic type to be used within the CasADi functions of the object
            self.sym = 'SX';
            self.model_type = model_type;
            % User can now choose whether to input new data or load existing.
            %
            % Default numeric data is set via purpose template functions, we keep separate from class definitions.            
            % TODO Ask her.
            switch model_type
                case 'EVMultibodyVehicle'
                    self.data = ev_multibody_vehicle_data_template();
                case 'FSAEVehicle'
                    self.data = fsae_multibody_vehicle_data_template(track, alfarange);
                case 'DoubleTrackVehicle'
                    self.data = double_track_vehicle_data_template();
                case 'SingleTrackVehicle'
                    self.data = single_track_vehicle_data_template();
            end
            
            % Build variables and bounds 
            self.nx = 24;
            self.Xb = struct( ...
                'lb',[ 0;-1;-1;-1;-1;-1;-1;-1;-1;-1; 0; 0; 0; 0;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1], ...
                'ub',[ 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1] ...
                );
            self.nu = 3;
            self.Ub = struct( ...
                'lb',[ 0;-1;-1], ...
                'ub',[ 1; 0; 1] ...
                );
            self.nuc = 0;
            self.nz = 15;
            self.Zb = struct( ...
                'lb',[-1; 0; 0;-1;-1;-1;-1;-1;-1;-1;-1; 0; 0; 0; 0], ...
                'ub',[ 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1] ...
                );
            self.nzc = 0;
            self.nwg = 0;
            self.neq = 0;       
            % Build car
            [~,self.dyn_fun,self.out_fun] = self.dynamic_model();
            
        end
        function [sys_fun,dyn_fun,out_fun] = dynamic_model(self)
            
            % Implement the Articulated Body Algorithm to simulate the dynamics of the multibody vehicle.
            %
            % This function creates CasADi functions that contain the main functionality of the model:
            %
            % [dt_x,v,hgt] = sys_fun(x,u,track_data)
            % dt_x = dyn_fun(x,u,track_data,v)
            % v = out_fun(x,u,track_data)
            %
            % State variables: x_unscaled = [V_gb;dt_trim;omega;q_gb;trim]
            % Input variables: u_unscaled = [Ta;Tb;delta]
            % Sparsifying variables: v_unscaled = [Fx;Fy;Fz]
            % Time-varying parameters: track_data = [d_gs;reshape(R_gs,9,1)]
            % Additional info required for animation: hgt = {g_gb,g_gh_11,g_gh_12,g_gh_21,g_gh_22}
            
            import casadi.*
            
            %% Preliminary operations
            
            % Extract the relevant data
            q_bc{1,1} = self.data.q_bc_FL(1:6);
            q_bc{1,2} = self.data.q_bc_FL(1:6).*[ 1;-1; 1;-1; 1;-1];
            q_bc{2,1} = self.data.q_bc_RL(1:6);
            q_bc{2,2} = self.data.q_bc_RL(1:6).*[ 1;-1; 1;-1; 1;-1];
            a1 = self.data.a1;
            a2 = self.data.a2;
            m_w = self.data.m_w;
            i_xx_w = self.data.i_xx_w;
            i_yy_w = self.data.i_yy_w;
            m_b = self.data.m_b;
            i_xx_b = self.data.i_xx_b;
            i_yy_b = self.data.i_yy_b;
            i_zz_b = self.data.i_zz_b;
            i_xz_b = self.data.i_xz_b;
            tire_radius = self.data.tire_radius;
            brake_balance = self.data.brake_balance;
            h0 = self.data.h0;
            g = self.data.g;
            x_scale = self.data.X_scale;
            u_scale = self.data.U_scale;
            v_scale = self.data.Z_scale(4:15);

            % Define CasADi symbolic variables
            x = casadi.(self.sym).sym('x',self.nx,1);
            u = casadi.(self.sym).sym('u',self.nu,1);
            v = casadi.(self.sym).sym('v',12,1);
            track_data = casadi.(self.sym).sym('track_data',12,1);
            
            % Apply scaling factors
            x_unscaled = x_scale.*x;
            u_unscaled = u_scale.*u;
            v_unscaled = v_scale.*v;
            
            % Unpack states
            V_gb = x_unscaled(1:6);
            dt_trim{1,1} = x_unscaled(7);
            dt_trim{1,2} = x_unscaled(8);
            dt_trim{2,1} = x_unscaled(9);
            dt_trim{2,2} = x_unscaled(10);
            omega{1,1} = x_unscaled(11);
            omega{1,2} = x_unscaled(12);
            omega{2,1} = x_unscaled(13);
            omega{2,2} = x_unscaled(14);
            q_gb = x_unscaled(15:20);
            trim{1,1} = x_unscaled(21);
            trim{1,2} = x_unscaled(22);
            trim{2,1} = x_unscaled(23);
            trim{2,2} = x_unscaled(24);
            
            % Unpack inputs
            Ta = u_unscaled(1);
            Tb = u_unscaled(2);
            delta = u_unscaled(3);
            
            % Unpack sparsifying variables
            Fx{1,1} = v_unscaled(1);
            Fx{1,2} = v_unscaled(2);
            Fx{2,1} = v_unscaled(3);
            Fx{2,2} = v_unscaled(4);
            Fy{1,1} = v_unscaled(5);
            Fy{1,2} = v_unscaled(6);
            Fy{2,1} = v_unscaled(7);
            Fy{2,2} = v_unscaled(8);
            Fz{1,1} = v_unscaled(9);
            Fz{1,2} = v_unscaled(10);
            Fz{2,1} = v_unscaled(11);
            Fz{2,2} = v_unscaled(12);
            
            % Extract track data
            d_gs = track_data(1:3);
            R_gs = reshape(track_data(4:12),3,3);
            
            % Extract some convenience quantities from the configuration vector
            d_gb = q_gb(1:3);
            R_gb = rotZYX(q_gb(4:6));
            
            % Some symplifying assumptions
            dt_delta = 0;
            dt_dt_delta = 0;
            
            % Recall: cell indices {i:1|2,j:1|2} refer to the {axle:F|R,side:L|R} wheel.
            % Recall: unscaled is dimensional, scaled is adimensional (i.e. normalized, typically in [-1,1]).
            %
            % Now we are ready to apply the Articulated Body Algorithm.
            %
            % To start, recall the general notation for frame {B} w.r.t. {A}:
            % d_ab = origin position
            % E_ab = Euler angles
            % q_ab = 6D rigid-body coordinates (origin position + Euler angles)
            % R_ab = 3x3 rotation matrix
            % g_ab = 4x4 homogeneous transformation matrix
            %
            % Some remarks about the screw notation used in this function:
            % V_b = rigid-body velocity of {B} w.r.t. {G} expressed in {B}
            % V_h = rigid-body velocity of {H} w.r.t. {G} expressed in {H}
            % V_w = rigid-body velocity of {W} w.r.t. {G} expressed in {H}
            % M, B, W = same frame and body as the corresponding V
            %
            % For more detailed reference, see https://www.mdpi.com/2411-9660/7/3/65.
            
            %% Evaluate suspensions
            
            % Fetch the CasADi functions for suspension analysis
            [susp_fun,tau_fun] = self.suspension_analysis();
            
            d_bh = cell(2,2);
            R_bh = cell(2,2);
            J_bh_trim = cell(2,2);
            J_bh_delta = cell(2,2);
            dt_J_bh_trim = cell(2,2);
            dt_J_bh_delta = cell(2,2);
            Ad_g_bh_star = cell(2,2);
            for i = 1:2
                for j = 1:2
                    
                    % Compute the kinematic quantities related to suspensions
                    [q_ch,J_bh_trim{i,j},J_bh_delta{i,j},dt_J_bh_trim{i,j},dt_J_bh_delta{i,j}] = susp_fun{i,j}(trim{i,j},delta,dt_trim{i,j},dt_delta);
                    d_bh{i,j} = q_bc{i,j}(1:3)+q_ch(1:3);
                    R_bh{i,j} = rotZYX(q_ch(4:6));
                    Ad_g_bh_star{i,j} = adjointStar(rtm(R_bh{i,j},d_bh{i,j}));
                    
                end
            end
            
            % Define the Jacobian for the hub bearing
            J_hw = [0;0;0;0;1;0];
            
            tau = cell(2,2);
            for i = 1:2
                
                % Compute the generalized suspension force
                tau{i,1} = tau_fun{i}(trim{i,1},dt_trim{i,1},trim{i,2});
                tau{i,2} = tau_fun{i}(trim{i,2},dt_trim{i,2},trim{i,1});
                
            end
            
            %% Sweep I
            
            % Get the twist velocity of the chassis, already available from the state vector (here we'll use the short notation)
            V_b = V_gb;
            
            V_bh = cell(2,2);
            V_h = cell(2,2);
            V_hw = cell(2,2);
            V_w = cell(2,2);
            for i = 1:2
                for j = 1:2
                    
                    % Compute the twist velocity of knuckles
                    V_bh{i,j} = J_bh_trim{i,j}*dt_trim{i,j}+J_bh_delta{i,j}*dt_delta;
                    V_h{i,j} = Ad_g_bh_star{i,j}.'*V_b+V_bh{i,j};
                    
                    % Compute the twist velocity of rims
                    V_hw{i,j} = J_hw*omega{i,j};
                    V_w{i,j} = V_h{i,j}+V_hw{i,j};
                    
                end
            end
            
            %% Compute the external wrenches
            
            % Allocate the input torque to each wheel
            T_drive{1,1} = 0;
            T_drive{1,2} = 0;
            T_drive{2,1} = 0.5*Ta;
            T_drive{2,2} = 0.5*Ta;
            T_brake{1,1} = 0.5*brake_balance*Tb*tanh(1*omega{1,1});
            T_brake{1,2} = 0.5*brake_balance*Tb*tanh(1*omega{1,2});
            T_brake{2,1} = 0.5*(1-brake_balance)*Tb*tanh(1*omega{2,1});
            T_brake{2,2} = 0.5*(1-brake_balance)*Tb*tanh(1*omega{2,2});
            
            % Fetch the CasADi function for tire force evaluation
            tire_fun = self.tire_model();
            
            % Preliminarily express chassis w.r.t. track
            R_sb = R_gs.'*R_gb;
            
            W_h_ext = cell(2,2);
            W_w_ext = cell(2,2);
            R_gh = cell(2,2);
            R_sh = cell(2,2);
            d_hn = cell(2,2);
            R_hn = cell(2,2);
            Fx_out = cell(2,2);
            Fy_out = cell(2,2);
            Fz_out = cell(2,2);
            for i = 1:2
                for j = 1:2
                    
                    % Express wheel w.r.t. road
                    R_gh{i,j} = R_gb*R_bh{i,j};
                    R_sh{i,j} = R_sb*R_bh{i,j};
                    
                    % Compute the wheel angles w.r.t. road
                    gamma = atan2(R_sh{i,j}(3,2),sqrt(1-R_sh{i,j}(3,2)^2)); % NOTE The camber angle.
                    sigma = atan2(-R_sh{i,j}(3,1),R_sh{i,j}(3,3)); % NOTE The pitch angle.
                    
                    % Compute the height of wheel center above the road surface
                    h = R_gs(:,3).'*(R_gb*d_bh{i,j}+d_gb-d_gs);
                    
                    % Express ground (aligned to tire) w.r.t. wheel
                    d_hn{i,j} = h*[sin(sigma)/cos(gamma);0;-cos(sigma)/cos(gamma)];
                    R_hn{i,j} = rotY(-sigma)*rotX(-gamma);
                    
                    % Compute the tire slips
                    V = V_w{i,j}(1:3)+cross(V_w{i,j}(4:6),d_hn{i,j}); % NOTE Velocity of tire contact point w.r.t. wheel.
                    V = R_hn{i,j}.'*V; % NOTE Velocity of tire contact point w.r.t. road.
                    Vx = R_hn{i,j}(:,1).'*V_h{i,j}(1:3);
                    den = abs(Vx)+eps*(Vx==0); % TODO Fix this.
                    kappa = V(1)/den;
                    alpha = -V(2)/den;
                    
                    % Compute the tire displacement for the penalty contact model
                    d = tire_radius{i}*cos(gamma)-h;
                    dt_d = -V(3);
                    
                    % Compute the components of the tire force (to be directly output for sparsification)
                    [Fx_out{i,j},Fy_out{i,j},Fz_out{i,j}] = tire_fun{i}(kappa,alpha,gamma,d,dt_d);
                    
                    % Compute the weight force on the wheel
                    Fg_w = -m_w{i}*g;
                    
                    % Compute the wrench on the rims (directly using the sparsifying variables)
                    W_w_ext{i,j} = adjointStar(rtm(R_hn{i,j},d_hn{i,j}))*[Fx{i,j};Fy{i,j};Fz{i,j};0;0;0] ... % NOTE The contribution of the tire force.
                        +[R_gh{i,j}.'*[0;0;Fg_w];zeros(3,1)] ... % NOTE The contribution of the weight force.
                        +J_hw*(T_drive{i,j}+T_brake{i,j}); % NOTE The contribution of the input torque.
                    
                    % Compute the wrench on the knuckles
                    W_h_ext{i,j} = -J_hw*T_brake{i,j}; % NOTE The reaction of the braking torque.
                    
                end
            end
            
            % Fetch the CasADi function for the aerodynamic force
            aero_fun = self.aerodynamic_model();
            
            % Compute the aerodynamic forces
            [Fd,Fl1,Fl2] = aero_fun(V_b(1));
            
            % Compute the weight force on the chassis
            Fg_b = -m_b*g;
            
            % Compute the wrench on the chassis
            W_b_ext = [-Fd;0;-Fl1-Fl2;0;a1*Fl1-a2*Fl2+h0*Fd;0] ... % NOTE The contribution of the aerodynamic force.
                +[R_gb.'*[0;0;Fg_b];zeros(3,1)] ... % NOTE The contribution of the weight force.
                -Ad_g_bh_star{1,1}*J_hw*T_drive{1,1} ...
                -Ad_g_bh_star{1,2}*J_hw*T_drive{1,2} ...
                -Ad_g_bh_star{2,1}*J_hw*T_drive{2,1} ...
                -Ad_g_bh_star{2,2}*J_hw*T_drive{2,2}; % NOTE The reactions of the driving torques.
            
            %% Sweep II
            
            M_w_hat = cell(2,2);
            B_w_hat = cell(2,2);
            M_h_hat = cell(2,2);
            B_h_hat = cell(2,2);
            K = cell(2,2);
            I = cell(2,2);
            M = cell(2,2);
            U = cell(2,2);
            T = cell(2,2);
            B = cell(2,2);
            for i = 1:2
                for j = 1:2
                    
                    % Compute the articulated body inertia of rims
                    M_w_hat{i,j} = diag([m_w{i},m_w{i},m_w{i},i_xx_w{i},i_yy_w{i},i_xx_w{i}]);
                    
                    % Compute the articulated body bias of rims
                    B_w_hat{i,j} = -W_w_ext{i,j}+adStar(V_w{i,j})*(M_w_hat{i,j}*V_w{i,j});
                    
                    % Compute the articulated body inertia of knuckles
                    M_h_hat{i,j} = diag([m_w{i},m_w{i},m_w{i},i_xx_w{i},0,i_xx_w{i}]);
                    
                    % Compute the articulated body bias of knuckles
                    B_h_hat{i,j} = -W_h_ext{i,j}+B_w_hat{i,j}-J_hw*B_w_hat{i,j}(5)-M_h_hat{i,j}*ad(V_hw{i,j})*V_w{i,j};
                    
                    % Compute some common sub-expressions (their names are arbitrary, and have no special meaning)
                    K{i,j} = M_h_hat{i,j}*J_bh_trim{i,j};
                    I{i,j} = K{i,j}.'*J_bh_trim{i,j};
                    M{i,j} = M_h_hat{i,j}-K{i,j}*(K{i,j}.'/I{i,j});
                    U{i,j} = dt_J_bh_trim{i,j}*dt_trim{i,j}+dt_J_bh_delta{i,j}*dt_delta+J_bh_delta{i,j}*dt_dt_delta-ad(V_bh{i,j})*V_h{i,j};
                    T{i,j} = tau{i,j}-J_bh_trim{i,j}.'*B_h_hat{i,j};
                    B{i,j} = B_h_hat{i,j}+M{i,j}*U{i,j}+T{i,j}/I{i,j}*K{i,j};
                    
                end
            end
            
            % Compute the articulated body inertia of the chassis
            M_b = diag([m_b,m_b,m_b,i_xx_b,i_yy_b,i_zz_b]);
            M_b(4,6) = i_xz_b;
            M_b(6,4) = i_xz_b;
            M_b_hat = M_b ...
                +Ad_g_bh_star{1,1}*M{1,1}*Ad_g_bh_star{1,1}.' ...
                +Ad_g_bh_star{1,2}*M{1,2}*Ad_g_bh_star{1,2}.' ...
                +Ad_g_bh_star{2,1}*M{2,1}*Ad_g_bh_star{2,1}.' ...
                +Ad_g_bh_star{2,2}*M{2,2}*Ad_g_bh_star{2,2}.';
            
            % Compute the articulated body inertia of the chassis
            B_b_hat = -W_b_ext+adStar(V_b)*(M_b*V_b) ...
                +Ad_g_bh_star{1,1}*B{1,1} ...
                +Ad_g_bh_star{1,2}*B{1,2} ...
                +Ad_g_bh_star{2,1}*B{2,1} ...
                +Ad_g_bh_star{2,2}*B{2,2};
            
            %% Sweep III
            
            % Compute the twist acceleration of the chassis
            dt_V_b = -M_b_hat\B_b_hat;
            
            Y = cell(2,2);
            dt_dt_trim = cell(2,2);
            dt_V_h = cell(2,2);
            dt_omega_ij = cell(2,2);
            for i = 1:2
                for j = 1:2
                    
                    % Compute the suspension travel acceleration
                    Y{i,j} = Ad_g_bh_star{i,j}.'*dt_V_b+U{i,j};
                    dt_dt_trim{i,j} = (T{i,j}-K{i,j}.'*Y{i,j})/I{i,j};
                    
                    % Compute the twist acceleration of knuckles
                    dt_V_h{i,j} = Y{i,j}+J_bh_trim{i,j}*dt_dt_trim{i,j};
                    
                    % Compute the wheel acceleration
                    dt_omega_ij{i,j} = -B_w_hat{i,j}(5)/M_w_hat{i,j}(5,5)-dt_V_h{i,j}(5);
                    
                end
            end
            
            %% Compute derivative of linear position and Euler angles w.r.t. world
            
            dt_q_gb = floatingBaseBodyJacInvZYX(q_gb)*V_gb;
            
            %% Save CasADi functions
            
            dt_V_gb = dt_V_b;
            
            % Pack state derivatives
            dt_x_unscaled = [
                dt_V_gb;
                dt_dt_trim{1,1};
                dt_dt_trim{1,2};
                dt_dt_trim{2,1};
                dt_dt_trim{2,2};
                dt_omega_ij{1,1};
                dt_omega_ij{1,2};
                dt_omega_ij{2,1};
                dt_omega_ij{2,2};
                dt_q_gb;
                dt_trim{1,1};
                dt_trim{1,2};
                dt_trim{2,1};
                dt_trim{2,2};
                ];
            
            % Pack sparsifying variables
            v_out_unscaled = [
                Fx_out{1,1};
                Fx_out{1,2};
                Fx_out{2,1};
                Fx_out{2,2};
                Fy_out{1,1};
                Fy_out{1,2};
                Fy_out{2,1};
                Fy_out{2,2};
                Fz_out{1,1};
                Fz_out{1,2};
                Fz_out{2,1};
                Fz_out{2,2};
                ];
            
            % Reapply scaling factors
            dt_x = x_scale.\dt_x_unscaled;
            v_out = v_scale.\v_out_unscaled;
            
            % Prepare additional outputs
            g_gb = rtm(R_gb,d_gb);
            g_gh_11 = rtm(R_gh{1,1},R_gb*d_bh{1,1}+d_gb);
            g_gh_12 = rtm(R_gh{1,2},R_gb*d_bh{1,2}+d_gb);
            g_gh_21 = rtm(R_gh{2,1},R_gb*d_bh{2,1}+d_gb);
            g_gh_22 = rtm(R_gh{2,2},R_gb*d_bh{2,2}+d_gb);
            
            % Save CasADi functions
            out_fun = casadi.Function('out',{x,u,track_data},{v_out},struct('cse',true));
            dyn_fun = casadi.Function('dyn',{x,u,track_data,v},{dt_x},struct('cse',true));
            dt_x_out = dyn_fun(x,u,track_data,v_out);
            sys_fun = casadi.Function('dyn',{x,u,track_data},{dt_x_out,v_out,g_gb,g_gh_11,g_gh_12,g_gh_21,g_gh_22},struct('cse',true));
            
        end
        function tire_fun = tire_model(self)
            
            % Define the tire model (tangential friction + normal contact) to be used in multibody vehicle simulation.
            

            
            import casadi.*
            
            % Extract the relevant data
            tir = self.data.tir;
            
            % Define CasADi symbolic variables
            kappa = casadi.(self.sym).sym('kappa',1,1);
            alpha = casadi.(self.sym).sym('alpha',1,1);
            gamma = casadi.(self.sym).sym('gamma',1,1);
            d = casadi.(self.sym).sym('d',1,1);
            dt_d = casadi.(self.sym).sym('dt_d',1,1);
            
            tire_fun = cell(2,1);
            for i = 1:2
                
                % Get the relevant tire properties
                tire_k = self.data.tire_k{i};
                tire_c = self.data.tire_c{i};
                
                % Compute the vertical load
                Fz = (tire_k*d+tire_c*dt_d)*if_else_smooth(d,0,1,0,'C',1e3);
                
                % Fetch the Magic Formula CasADi function
                [~,~,~,~,~,~,~,~,mf_fun] = MF_Tyre(tir{i});
                %mf_fun = magic_formula(tir{i},'sym',self.sym);
                
                % Evaluate the Magic Formula
                [Fx,Fy] = mf_fun(kappa,alpha,Fz,gamma);
                
                % Save CasADi function for the front wheels
                tire_fun{i} = casadi.Function('tire',{kappa,alpha,gamma,d,dt_d},{Fx,Fy,Fz},struct('cse',true));
                
            end
            
        end
        function aero_fun = aerodynamic_model(self)
            
            % Define the aerodynamic model to be used in multibody vehicle simulation.
            
            import casadi.*
            
            % Extract the relevant data
            S = self.data.aero_S;
            rho = self.data.aero_rho;
            Cx = self.data.aero_Cx;
            Cz = self.data.aero_Cz;
            balance = self.data.aero_balance;

            % Define CasADi symbolic variables
            u = casadi.(self.sym).sym('u',1,1);
            
            % Compute the drag force
            Fd = 0.5*S*rho*Cx*u^2*(u>0);
            
            % Compute the downforce on the front axle
            Fl1 = balance*0.5*S*rho*Cz*u^2;
            
            % Compute the downforce on the rear axle
            Fl2 = (1-balance)*0.5*S*rho*Cz*u^2;
            
            % Save CasADi function
            aero_fun = casadi.Function('tire',{u},{Fd,Fl1,Fl2},struct('cse',true));
            
        end
        function [susp_fun,tau_fun] = suspension_analysis(self,options)
            
            % Solve the kinematics of double-wishbone (or generic 5-link) suspensions.
            
            arguments
                
                % This object
                self
                
                % Number of samples over the suspension travel range
                options.N_trim = 15;
                
                % Number of samples over the steering wheel range
                options.N_delta = 15;
                
                % Order of the regression polynomials
                options.deg = 3;

                
            end
            
            import casadi.*
            
            %% Preliminary operations
            
            % Extract the relevant data
            P_c_susp{1} = self.data.P_c_susp_FL;
            P_c_susp{2} = self.data.P_c_susp_RL;
            P_h_susp{1} = self.data.P_h_susp_FL;
            P_h_susp{2} = self.data.P_h_susp_RL;
            P_c_spring{1} = self.data.P_c_spring_FL;
            P_c_spring{2} = self.data.P_c_spring_RL;
            P_r_spring{1} = self.data.P_r_spring_FL;
            P_r_spring{2} = self.data.P_r_spring_RL;
            P_r_push{1} = self.data.P_r_push_FL;
            P_r_push{2} = self.data.P_r_push_RL;
            P_h_push{1} = self.data.P_h_push_FL;
            P_h_push{2} = self.data.P_h_push_RL;
            d_ch{1} = self.data.q_ch_FL(1:3);
            d_ch{2} = self.data.q_ch_RL(1:3);
            E_ch{1} = self.data.q_ch_FL(4:6);
            E_ch{2} = self.data.q_ch_RL(4:6);
            d_cr{1} = self.data.q_cr_FL(1:3);
            d_cr{2} = self.data.q_cr_RL(1:3);
            E_cr{1} = self.data.q_cr_FL(4:6);
            E_cr{2} = self.data.q_cr_RL(4:6);
            steer_dir = self.data.steer_dir;
            steer_ratio = self.data.steer_ratio;
            steer_max = self.data.steer_max;
            k_s = self.data.k_s;
            c_d = self.data.c_d;
            k_a = self.data.k_a;
            coil_length = self.data.coil_length;
            trim_lim = self.data.trim_lim;
            bump_order = self.data.bump_order;
            
            % Unpack the options
            N_trim = options.N_trim;
            N_delta = options.N_delta;
            deg = options.deg;
            
            % IPOPT settings
            opts.ipopt.print_level = 0;
            opts.print_time = 0;
            
            % Compute the rod lengths (they are computed at the given reference configuration)
            l_susp = cell(2,1);
            l_push = cell(2,1);
            l_spring = cell(2,1);
            for i = 1:2
                P_c_susp_c = P_c_susp{i};
                P_h_susp_c = rotZYX(E_ch{i})*P_h_susp{i}+d_ch{i};
                P_r_push_c = rotZYX(E_cr{i})*P_r_push{i}+d_cr{i};
                P_h_push_c = rotZYX(E_ch{i})*P_h_push{i}+d_ch{i};
                P_c_spring_c = P_c_spring{i};
                P_r_spring_c = rotZYX(E_cr{i})*P_r_spring{i}+d_cr{i};
                l_susp{i} = sqrt(sum((P_c_susp_c-P_h_susp_c).^2));
                l_push{i} = sqrt(sum((P_r_push_c-P_h_push_c).^2));
                l_spring{i} = sqrt(sum((P_c_spring_c-P_r_spring_c).^2)); % NOTE The spring length is variable, this is its value at reference.
            end
            
            % Define CasADi symbolic variables
            trim_sym = casadi.(self.sym).sym('trim',1,1);
            delta_sym = casadi.(self.sym).sym('delta',1,1);
            dt_trim_sym = casadi.(self.sym).sym('dt_trim',1,1);
            dt_delta_sym = casadi.(self.sym).sym('dt_delta',1,1);
            d_ch_sym = casadi.(self.sym).sym('d_ch',3,1);
            E_ch_sym = casadi.(self.sym).sym('E_ch',3,1);
            beta_sym = casadi.(self.sym).sym('beta',1,1);
            trim_other_sym = casadi.(self.sym).sym('trim_other',1,1);
            
            %% Fitting procedure
            
            % Delta is the steering wheel angle (independently controleld by the driver).
            % The trim is defined as the deviation of the spring length from its value at reference configuration.
            %
            % The fitting procedure is now performed for the FL and RL wheels.
            % The results are then mirrored to the FR and RR wheels.
            
            % Create monomial bases for the polynomials and their derivatives
            trim_mon = trim_sym.^(0:deg);
            d_trim_mon = (0:deg).*[0,trim_sym.^(0:deg-1)];
            dd_trim_mon = [0,(0:deg-1)].*(0:deg).*[0,0,trim_sym.^(0:deg-2)];
            delta_mon = (delta_sym.^(0:deg)).';
            d_delta_mon = ((0:deg).*[0,delta_sym.^(0:deg-1)]).';
            dd_delta_mon = ([0,(0:deg-1)].*(0:deg).*[0,0,delta_sym.^(0:deg-2)]).';
            
            % Sample the range of the independent variables (if range is symmetrical, midpoint is reference)
            trim_span{1} = linspace(trim_lim{1}(1),trim_lim{1}(2),N_trim);
            trim_span{2} = linspace(trim_lim{2}(1),trim_lim{2}(2),N_trim);
            delta_span = linspace(-steer_max,steer_max,N_delta);
            
            % Create the regression matrices
            trim_mat{1} = bsxfun(@power,trim_span{1}.',0:deg);
            trim_mat{2} = bsxfun(@power,trim_span{2}.',0:deg);
            delta_mat = bsxfun(@power,delta_span.',0:deg);
            
            susp_fun = cell(2,2);
            for i = 1:2
                
                %% Sample data for fitting
                
                % To relate the susp variables [d_ch;E_ch] and to trim and delta, we fit the data with 2D regression polynomials.
                %
                % The value of [d_ch;E_ch] first is computed for each combination of trim and delta.
                % These values are found upon evaluating the suspension kinematics.
                
                % Save rack endpoint position at reference configuration
                P_c_susp_5 = P_c_susp{i}(:,5);
                
                q_ch_data = nan(N_trim,N_delta,6);
                for ii = 1:N_trim
                    
                    % Find the value of the rocker angle corresponding to current trim
                    R_cr_sym = rotZYX(E_cr{i})*rotZ(beta_sym);
                    P_r_spring_c_sym = R_cr_sym*P_r_spring{i}+d_cr{i};
                    constr_beta = sqrt(sum((P_c_spring{i}-P_r_spring_c_sym).^2))-(l_spring{i}+trim_span{i}(ii));
                    problem_beta = struct('x',beta_sym,'f',0,'g',constr_beta);
                    solver_beta = nlpsol('solver','ipopt',problem_beta,opts);
                    beta_tmp = solver_beta('x0',0,'ubg',0,'lbg',0,'lbx',-pi/2,'ubx',+pi/2);
                    beta_tmp = full(beta_tmp.x);
                    
                    % Locate the push-rod mounting point on the rocker
                    R_cr_tmp = rotZYX(E_cr{i})*rotZ(beta_tmp);
                    P_r_push_c_tmp = R_cr_tmp*P_r_push{i}+d_cr{i};
                    
                    for jj = 1:N_delta
                        
                        % Move rack endpoint according to current delta
                        P_c_susp{i}(:,5) = P_c_susp_5+steer_dir{i}*steer_ratio{i}*delta_span(jj);
                        
                        % Write the equations of constraint and solve the suspension kinematics to find the values for the current sample
                        R_ch_sym = rotZYX(E_ch_sym);
                        P_h_susp_c_sym = R_ch_sym*P_h_susp{i}+d_ch_sym;
                        P_h_push_c_sym = R_ch_sym*P_h_push{i}+d_ch_sym;
                        constr_q_ch = [
                            sqrt(sum((P_c_susp{i}-P_h_susp_c_sym).^2))-l_susp{i}, ...
                            sqrt(sum((P_r_push_c_tmp-P_h_push_c_sym).^2))-l_push{i}, ...
                            ].'; % NOTE Last component is relative to the push-rod.
                        %                         constr_susp = (sqrt(sum(([P_c_susp{i},P_r_push_c_tmp]-[P_h_susp_c_sym,P_h_push_c_sym]).^2))-[l_susp{i},l_push{i}]).';
                        problem_q_ch = struct('x',[d_ch_sym;E_ch_sym],'f',0,'g',constr_q_ch);
                        solver_q_ch = nlpsol('solver','ipopt',problem_q_ch,opts);
                        q_ch_tmp = solver_q_ch('x0',[d_ch{i};E_ch{i}],'ubg',0,'lbg',0,'lbx',[d_ch{i}-1;E_ch{i}-pi],'ubx',[d_ch{i}+1;E_ch{i}+pi]);
                        q_ch_tmp = full(q_ch_tmp.x);
                        q_ch_data(ii,jj,:) = q_ch_tmp;
                        
                    end
                    
                end
                
                %% Fit data
                
                % The sampled data (for each component) are fitten with regression polynomials.
                % Polynomials are evaluated as [1 trim ... trim^d]*coeff*[1 delta ... delta^d]^T.
                
                % Compute matrix pseudoinverse
                trim_pinv = pinv(trim_mat{i});
                delta_pinv = pinv(delta_mat);
                
                q_ch_sym = casadi.(self.sym)(6,1);
                dtrim_q_ch_sym = casadi.(self.sym)(6,1);
                ddelta_q_ch_sym = casadi.(self.sym)(6,1);
                dtrimdtrim_q_ch_sym = casadi.(self.sym)(6,1);
                dtrimddelta_q_ch_sym = casadi.(self.sym)(6,1);
                ddeltaddelta_q_ch_sym = casadi.(self.sym)(6,1);
                for idx = 1:6
                    
                    % Compute the coefficients of the regression polynomials for the susp coordinates
                    coeff = trim_pinv*q_ch_data(:,:,idx)*delta_pinv.';
                    if i==2
                        coeff(:,2:end) = 0; % NOTE Only for a FWD vehicle.
                    end
                    
                    % Create CasADi expressions for the susp coordinates and their derivatives (componentwise)
                    q_ch_sym(idx) = trim_mon*coeff*delta_mon;
                    dtrim_q_ch_sym(idx) = d_trim_mon*coeff*delta_mon;
                    ddelta_q_ch_sym(idx) = trim_mon*coeff*d_delta_mon;
                    dtrimdtrim_q_ch_sym(idx) = dd_trim_mon*coeff*delta_mon;
                    dtrimddelta_q_ch_sym(idx) = d_trim_mon*coeff*d_delta_mon;
                    ddeltaddelta_q_ch_sym(idx) = trim_mon*coeff*dd_delta_mon;
                    
                end
                
                %% Compute Jacobians
                
                % Knowing the susp variables, we can compute the suspension Jacobians.
                
                % Compute the susp Jacobians
                J_ch_sym = floatingBaseBodyJacZYX(q_ch_sym);
                J_ch_trim_sym = J_ch_sym*dtrim_q_ch_sym;
                J_ch_delta_sym = J_ch_sym*ddelta_q_ch_sym;
                
                % Compute the derivative of the susp Jacobians
                dt_q_ch_sym = dtrim_q_ch_sym*dt_trim_sym+ddelta_q_ch_sym*dt_delta_sym;
                dt_J_ch_sym = dt_floatingBaseBodyJacZYX(J_ch_sym,dt_q_ch_sym);
                dt_J_ch_trim_sym = dt_J_ch_sym*dtrim_q_ch_sym+J_ch_sym*(dtrimdtrim_q_ch_sym*dt_trim_sym+dtrimddelta_q_ch_sym*dt_delta_sym);
                dt_J_ch_delta_sym = dt_J_ch_sym*ddelta_q_ch_sym+J_ch_sym*(dtrimddelta_q_ch_sym*dt_trim_sym+ddeltaddelta_q_ch_sym*dt_delta_sym);
                
                %% Save CasADi function
                
                % Save CasADi function for the left wheel
                susp_fun{i,1} = casadi.Function( ...
                    'susp', ...
                    {trim_sym,delta_sym,dt_trim_sym,dt_delta_sym}, ...
                    {q_ch_sym,J_ch_trim_sym,J_ch_delta_sym,dt_J_ch_trim_sym,dt_J_ch_delta_sym}, ...
                    struct('cse',true) ...
                    );
                
                % Mirror results for the right wheel
                [q_ch_sym,J_ch_trim_sym,J_ch_delta_sym,dt_J_ch_trim_sym,dt_J_ch_delta_sym] = susp_fun{i,1}(trim_sym,-delta_sym,dt_trim_sym,-dt_delta_sym);
                q_ch_sym = q_ch_sym.*[ 1;-1; 1;-1; 1;-1]; % NOTE These mirroring coefficients depend on the Euler sequence we are using.
                J_ch_trim_sym = J_ch_trim_sym.*[ 1;-1; 1;-1; 1;-1]; % NOTE These coefficients may be different from the coeff above.
                J_ch_delta_sym = J_ch_delta_sym.*[ 1;-1; 1;-1; 1;-1];
                dt_J_ch_trim_sym = dt_J_ch_trim_sym.*[ 1;-1; 1;-1; 1;-1];
                dt_J_ch_delta_sym = dt_J_ch_delta_sym.*[ 1;-1; 1;-1; 1;-1];
                
                % Save CasADi function for the right wheel
                susp_fun{i,2} = casadi.Function( ...
                    'susp', ...
                    {trim_sym,delta_sym,dt_trim_sym,dt_delta_sym}, ...
                    {q_ch_sym,J_ch_trim_sym,J_ch_delta_sym,dt_J_ch_trim_sym,dt_J_ch_delta_sym}, ...
                    struct('cse',true) ...
                    );
                
            end
            
            %% Compute the suspension force
            
            % The generalized suspension force relative to the trim is just the force along the spring.
            
            tau_fun = cell(2,1);
            for i = 1:2
                
                % Compute the generalized suspension force
                tau_sym = ( ...
                    -k_s{i}*(trim_sym+l_spring{i}-coil_length{i}) ... % NOTE The contribution of the spring.
                    -c_d{i}*dt_trim_sym ... % NOTE The contribution of the damper.
                    -k_a{i}*(trim_sym-trim_other_sym) ... % NOTE The contribution of the anti-roll bar.
                    -((2*trim_sym-trim_lim{i}(2)-trim_lim{i}(1))/abs(trim_lim{i}(2)-trim_lim{i}(1)))^bump_order{i} ... % NOTE The contribution of the hard-stops.
                    );
                
                % Save CasADi function
                tau_fun{i} = casadi.Function('tau',{trim_sym,dt_trim_sym,trim_other_sym},{tau_sym},struct('cse',true));
                
            end
            
        end
        function [x,v] = simulate(self,t_start,t_stop,h,x_init,u,track_data)
            
            % Run offline simulation.
            %
            % This function returns the (scaled) state trajectory and the value of the tire forces. We may add some additional outputs later.
            
            arguments
                self
                t_start
                t_stop
                h
                x_init
                u
                track_data = nan
            end
            import casadi.*
            % This function uses the default CasADi symbolic type for this vehicle
            
            % Fetch the system functions
            sys_fun = self.dynamic_model();
            
            % Discretize for off-line simulation
            x_sym = casadi.(self.sym).sym('x',self.nx,1);
            u_sym = casadi.(self.sym).sym('u',self.nu,1);
            track_data_sym = casadi.(self.sym).sym('track_data',12,1);
            h_sym = h;
            x_next_sym = RK4_step(x_sym,u_sym,@(x,u) sys_fun(x,u,track_data_sym),h_sym); % TODO Enable the user to choose the integrator. Here RK4 is used by default.
            [~,v_sym] = sys_fun(x_sym,u_sym,track_data_sym);
            next_step_fun = casadi.Function('next_step',{x_sym,u_sym,track_data_sym},{x_next_sym,v_sym},struct('cse',true));
            
            % Compile
%             opts = struct('mex',true,'main',false);
%             next_step_fun.generate('next_step_mex.c',opts);
%             mex next_step_mex.c -largeArrayDims
            
            % Define the time base
            t = t_start:h:t_stop;
            N = length(t)-1;
            
            % Define track data (if not provided). By default we run on a flat ground.
            if isnan(track_data)
                track_data = repmat([zeros(3,1);reshape(eye(3),9,1)],1,N);
            end
            
            % Simulate
            x = nan(self.nx,N+1);
            v = nan(12,N);
            x(:,1) = x_init;
            for k = 1:N
                [x_next,v_next] = next_step_fun(x(:,k),u(:,k),track_data(:,k));
                x(:,k+1) = full(x_next);
                v(:,k) = full(v_next);
            end
            
        end

    end
    

end