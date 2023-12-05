classdef track_vechicle_model < handle

    properties % visible properties from workspace
        
        track           % track where the vehicle is modeled on
        q0              % initial configuration of the vechicle w.r.t. track frame
        C               % daming coefficients array (d, theta, phi) 
        K               % stiffness coefficients array (d, theta, phi)
        xdot_sym        % symbolic function of state-space dynamics
        axle_wrench_sym % symbolic function of the dynamic wrench on axle
        axle_twist      % symbolic function of the axle twist w.r.t. ground
        xdot_subGraph
        axle_wrench_subGraph
        axle_twist_subGraph
        q               % positional states [alpha, ep, epsi, d, theta, phi]
        qdot            % positional states derivatives [alpha, ep, epsi, d, theta, phi]
        V               % axle absolute twist in axle frame
        w2              % wrench in axle frame {2} [Fx; Fy; Fz; Mx; My; Mz]
        iteration = 0   % integration iteration counter
        data            % store all extra data 
        fun             % store all extra function
    end
    
    properties (Hidden = true) % not visible properties
        
        B = eye(4);     % track-body transform matrix
        Ggb = eye(4);   % spatial-body transform matrix
        S = eye(4);     % spatial-track transform matrix
        X           % unit twist (in local frames)
        transformS  % transform handle for track frame
        transformB  % transform handle for body frame
        axes        % plot axes
        g0          % initial offset transform between local frames
        Blocal      % transform matrix of each joint in local coordinates
        bT          % twist of {S} w.r.t. {G} in {S} components
        M           % vehicle generalized inertia in principal baricentric frame {B}
        Maxle       % axle generalized inertia in baricentric frame {2}
        a1          %
        a2          %
        cx
        cz1
        cz2
        rho
        Sf
        GM          %
        k_phi1
        k_phi2
        k_phi
        q1
        q2
        qm
        h
        T_jac      % kinematic jacobian of track w.r.t. alpha (geometric twist w.r.t. to the nurbs variable alpha)
        T_jac_der  % derivative of track jacobian
        BS_solver = []; %solver for BS
        
    end
    
    properties (Access = private)
        
        
        
    end
    
    methods (Hidden = true) % hidden methods (hidden to user when promting hint)
        
        function J = track_to_body_jacobian(obj, q)
            %
            % {S} to {B} body jacobian expressed in {B}
            %
            
            J = BodyJac_localPOE(obj.g0, obj.X, q);
            
            
        end
        
        function J  = axle_to_body_jacobian(obj, q)
            %
            % {2} to {B} spatial jacobian expressed in {2}
            %
            
            J = spatialJac_localPOE(obj.g0(3:5), obj.X(:, 3:5), q(3:5));
            
        end
        
        function J = track_to_axle_jacobian(obj, q)
            %
            % {S} to {2} body jacobian expressed in {2}
            %
            % q here is [ep, ef, d, theta, phi]
            
            J = BodyJac_localPOE(obj.g0(1:2), obj.X(:, 1:2), q(1:2));
            
        end
        
        function [Jac_S, bT] = baseTwist(obj, alpha, alphadot)
            %
            % computes the base twist of track local frame wrt to the
            % ground
            %
            % AND the jacobian in axle frame {2}
            
            %
            % alpha_dot = 1./track.fun_normts(alpha).*s_dot
            
            % track induced geometric twist
            V_gs_s = obj.track.fun_Vgs_s(alpha);
            
            % track Jacobian {G} to {S} column in {S} frame
            Jac = V_gs_s.*obj.track.fun_normts(alpha);
            
            % track induced twist
            V_gs_s = V_gs_s.*obj.track.fun_normts(alpha).*alphadot;
            
            obj.bT = V_gs_s;
            
            % output the twist if needed
            if nargout > 0
                bT = V_gs_s;
                Jac_S = Jac;
            end
            
        end
        
        function Vdot_b = bodyDyn(obj, q, V_b, wp, wg_b, qd)
            %
            % INPUTS:
            %   q = configuration variables of just the actuated joints [d, theta, phi]
            %   qd = joint velocity variables of just the actuated joints [d_dot, theta_dot, phi_dot]
            %   Vgb = body-ground twist written in the axles frame {2}
            %   wp = pneumatic forces given as a wrench in {2}
            %
            % OUTPUT:
            %   Vdot_b = body twist derivative in {B}
            %
            
            % extract unitary twists
            Xu = obj.X;
            
            % spring and dampers forces
            tau = nan(3, 1, class(q));
            tau(3) = -q(3).*obj.K(3) - qd(3).*obj.C(3);  % rolling momentum
            tau(2) = -q(2).*obj.K(2) - qd(2).*obj.C(2);  % pitch momentum
            tau(1) = -q(1).*obj.K(1) - qd(1).*obj.C(1);  % shaking force
            
            % roll wrench
            w_roll = tau(3).*Xu(:, 5);
            
            % pitch wrench
            G = rigidInverse(obj.Blocal{5}); % G_5,4
            w_pitch = adjointStar(G)*Xu(:, 4).*tau(2);
            
            % shaking wrench
            G = G*rigidInverse(obj.Blocal{4}); % G_5,3
            w_d = adjointStar(G)*Xu(:, 3).*tau(1);
            
            % pneumatic wrench applied in {B}
            G = G*rigidInverse(obj.Blocal{3}); % G_5,2
            w_p = adjointStar(G)*wp;
            
            % total wrench applied to {B}
            w_b = w_roll + w_pitch + w_d + w_p + wg_b;
            
            % body equilibrium in {B} written in ODE form
            Vdot_b = obj.M\(w_b - adStar(V_b)*obj.M*V_b );
            
        end
        
        function eq = RNEAdyn(obj, alpha, alphadot, q, qd, qdd)
            %
            % not functional
            %
            
            % transform matrices update
            obj.FWKin_G(alpha, q);
            Xu = obj.X;
            
            % twist and acceleration forward propagation
            [Vbase, Jbase] = obj.baseTwist(alpha, alphadot); % base twist
            twist = nan(6, 6, class(q));  % twist init
            twistD = nan(6, 6, class(q)); % twist derivative initialization
            twist(:, 1) = Vbase; % base twist of the {S} frame
            twistD(:, 1) = [0;0;9.81;0;0;0]; % base acceleration (gravity)
            AdG = cell(1,5); % adjoint matrices init
            for i = 1:5
                AdG{i} = adjointInv(obj.Blocal{i});
                
                % collecting repeating terms for graph simplification
                AdTw = AdG{i}*twist(:, i);
                Xqd = Xu(:, i).*qd(i);
                
                twist(:, i+1) = AdTw + Xqd;
                twistD(:, i+1) = AdG{i}*twistD(:, i) + Xu(:, i).*qdd(i) + ad(AdTw)*Xqd;
                
            end
            
            % backwards propagation of wrenches
            % tau computation
            tau(1:2) = 0;
            tau(3) = -obj.K(1)*q(3) - obj.C(1)*qd(3);
            tau(4) = -obj.K(2)*q(4) - obj.C(2)*qd(4);
            tau(5) = -obj.K(3)*q(5) - obj.C(3)*qd(5);
            
            % inertia matrices
            II{1} = zeros(6, 6, class(q));
            II{2} = obj.Maxle; % axle inertia
            II{3} = II{1};
            II{4} = II{1};
            II{5} = obj.M; % body inertia
            
            eq = nan(1, 5, class(q)); % dynamic equilibrium equations in algebraic form
            w = nan(6,6, class(q)); % wrench init
            w(:, 6) = zeros(6, 1, class(q));
            AdG{i}{6} = adjointInv(eye(4)); % placeholder tranfrom for external wrench application
            
            for i = 5:1
                w(:, i) = AdG{i+1}.'*w(:, i+1) + II{i}*twistD(:, i+1) + adStar(twist(:, i+1))*II{i}*twist(:, i+1);
                eq(i) = Xu(:, i).'*w(:, i) - tau(i);
            end
            
        end
        
        function [qdd, twistAxle, wrenchAxle] = ABAdyn(obj, q, qd, wp) % tau is implicitly in the expression tau = fun(q, qd)
            %
            % [qdd, twistAxle, wrenchAxle] = ABAdyn(obj, Vbase, q, qd, wp)
            %
            % returns the virtual joint accelerations q_dotdotof the model
            % obtained with the Articulated Body Algorithm (ABA)
            %
            % states q = [alpha, ep, epsi, d, theta, phi]
            %        qdot = [alphadot, epdot, epsidot, ddot, thetadot, phidot]
            %
            % model input: wp: equivalent wrench applied to the axle frame
            %                 (third kinematic body)
            %
            % base twist is assumed null and the base acceleration is
            % assumed with an upward (along z) gravitational acceleration
            alpha = q(1);
            alphadot = qd(1);
            % transform matrices update
            obj.FWKin_G(alpha, q(2:6));
            % add track jacobian as an equivalent unitary twist
            Xu = [full(obj.T_jac(alpha)), obj.X];
            
            baseAcc = [0;0;9.81;0;0;0]; % simulating an upward gravitational acceleration
            
            % twist forward propagation
            twist = nan(6, 6, class(q));  % twist init
            twistD = nan(6, 6, class(q)); % twist derivative initialization
            AdGinv = cell(1, 6); % adjoint matrices init
            a = cell(6, 1);                % init relative twist ad
            
            % first joint (track curvilinear joint) is a special case
            twist(:, 1) = Xu(:, 1).*alphadot; % twist of the track frame
            AdGinv{1} = adjointInv(obj.S);
            %            a{1} = ad(Xu(:, 1).*alphadot);
            
            % propagation loop of twits
            for i = 2:6
                AdGinv{i} = adjointInv(obj.Blocal{i-1});
                
                twist(:, i) = AdGinv{i}*twist(:, i-1) + Xu(:, i).*qd(i); % twist of i-th body (i = 1, is the base)
                a{i} = ad(Xu(:, i).*qd(i));
                if i == 3
                    twistAxle = twist(:, i);
                end
            end
            
            % backward propagation of the articulated inertias and biases
            tau = nan(6, 1, class(q));
            tau(1:3) = 0;
            tau(4) = -obj.K(1)*(q(4)) - obj.C(1)*qd(4);  % vertical spring/damping
            tau(5) = -obj.K(2)*q(5) - obj.C(2)*qd(5);  % pitch torsional spring/damping
            tau(6) = -obj.K(3)*q(6) - obj.C(3)*qd(6);  % roll torsional spring/damping
            
            % inertia matrices
            if isa(q, 'casadi.SX'); zero = casadi.SX(6,6); else; zero = zeros(6,6); end
            II{1} = zero;
            II{2} = II{1};
            II{3} = obj.Maxle; % axle inertia
            II{4} = II{1};
            II{5} = II{1};
            II{6} = obj.M; % body inertia
            
            Mhat = cell(6, 1);             % init articulated inertia
            Mbar = cell(6, 1);             % init projected inertia
            bhat = nan(6, 6, class(q));    % init articulated bias
            Mscalar = nan(6, 1, class(q)); % init scalarazied non structural inertia (expression collection)
            
            for i = 6:-1:1
                if i == 6 % final body has no children
                    Mhat{i} = II{i};
                    Mscalar(i) = Xu(:, i).'*Mhat{i}*Xu(:, i) + 1e-3;
                    Mbar{i} = (  eye(6) - (Mhat{i}*Xu(:, i)*Xu(:, i).')./(Mscalar(i))  )*Mhat{i};
                    
                    bhat(:, i) = adStar(twist(:, i))*II{i}*twist(:, i);
                    g_axle_body = obj.Blocal{3}*obj.Blocal{4}*obj.Blocal{5};
                    rot_axle_body = g_axle_body(1:3,1:3);
                    rot_body_axle = rot_axle_body';
                    Faero_body = rot_body_axle*[-0.5*obj.cx*obj.rho*obj.Sf*twistAxle(1)^2,0,-0.5*(obj.cz1+obj.cz2)*obj.rho*obj.Sf*twistAxle(1)^2]';
                    Maero_body = rot_body_axle*[0,(obj.cz1*obj.a1-obj.cz2*obj.a2)*0.5*obj.rho*obj.Sf*twistAxle(1)^2,0]';
                    bhat(:, i) = bhat(:, i) - [Faero_body; Maero_body];
                    continue
                end
                
                j = i+1; % children body
                AdGstar = AdGinv{j}.';
                
                % current articulated inertia
                Mhat{i} = II{i} + AdGstar*Mbar{j}*AdGinv{j};
                Mscalar(i) = Xu(:, i).'*Mhat{i}*Xu(:, i) + 1e-3;
                
                % projected inertia (not required for the 1st body)
                if i ~= 1
                    Mbar{i} = (  eye(6) - (Mhat{i}*Xu(:, i)*Xu(:, i).')./(Mscalar(i))  )*Mhat{i};
                end
                
                % current body bias
                b_i = adStar(twist(:, i))*II{i}*twist(:, i);
                
                % current articulation bias
                V_ij = AdGinv{j}*twist(:, i); % current body twist in children frame
                bhat(:, i) = b_i + AdGstar*...
                    (...
                    bhat(:, j) - Mbar{j}*a{j}*V_ij + Mhat{j}*Xu(:, j).*(tau(j) - Xu(:, j).'*bhat(:, j))./Mscalar(j)...
                    ); % current articulated bias
                
                if i == 3 % axle body
                    bhat(:, i) = bhat(:, i) - wp;  % add pneumatic wrench
                end
            end
            
            % forward propagation of accelerations
            qdd = nan(6, 1, class(q));
            for i = 1:6
                if i == 1 % track joint transform
                    Vdot_parent = AdGinv{i}*baseAcc;
                    jac_acc = obj.T_jac_der(alpha).*alphadot.^2;
                    qdd(i) = -Xu(:, i).'*(Mhat{i}*(Vdot_parent + jac_acc) + bhat(:, i))./(Mscalar(i));
                    twistD(:, i) = Vdot_parent + Xu(:, i).*qdd(i) + jac_acc; % base velocity is always zero
                    continue
                end
                % current joint acceleration
                Vdot_parent = AdGinv{i}*twistD(:, i-1);
                V_p = a{i}*AdGinv{i}*twist(:, i-1);
                qdd(i) = (tau(i) - Xu(:, i).'*(Mhat{i}*(Vdot_parent - V_p) + bhat(:, i)))./(Mscalar(i));
                
                % current body twist derivative
                twistD(:, i) = Vdot_parent + Xu(:, i).*qdd(i) - V_p;
                
                % extract wrench on axles
                if i == 3
                    wrenchAxle = Mhat{i}*twistD(:, i) + bhat(:, i);
                end
            end
            
        end
        
        function [qdd, twistAxle, wrenchAxle] = ABAdyn_alpha_num(obj, g_gs_num, Jac_num, Jac_der_num, q, qd, wp) % tau is implicitly in the expression tau = fun(q, qd)
            %
            % [qdd, twistAxle, wrenchAxle] = ABAdyn(obj, Vbase, q, qd, wp)
            %
            % returns the virtual joint accelerations q_dotdotof the model
            % obtained with the Articulated Body Algorithm (ABA)
            %
            % states q = [alpha, ep, epsi, d, theta, phi]
            %        qdot = [alphadot, epdot, epsidot, ddot, thetadot, phidot]
            %
            % model input: wp: equivalent wrench applied to the axle frame
            %                 (third kinematic body)
            %
            % base twist is assumed null and the base acceleration is
            alphadot = qd(1);
            % transform matrices update
            obj.FWKin(q);
            obj.S = g_gs_num*obj.B;
            % add track jacobian as an equivalent unitary twist
            Xu = [Jac_num, obj.X];
            
            baseAcc = [0;0;9.81;0;0;0]; % simulating an upward gravitational acceleration
            
            % twist forward propagation
            twist = nan(6, 6, class(q));  % twist init
            twistD = nan(6, 6, class(q)); % twist derivative initialization
            AdGinv = cell(1, 6); % adjoint matrices init
            a = cell(6, 1);                % init relative twist ad
            
            % first joint (track curvilinear joint) is a special case
            twist(:, 1) = Xu(:, 1).*alphadot; % twist of the track frame
            AdGinv{1} = adjointInv(obj.S);
            %            a{1} = ad(Xu(:, 1).*alphadot);
            
            % propagation loop of twits
            for i = 2:6
                AdGinv{i} = adjointInv(obj.Blocal{i-1});
                
                twist(:, i) = AdGinv{i}*twist(:, i-1) + Xu(:, i).*qd(i); % twist of i-th body (i = 1, is the base)
                a{i} = ad(Xu(:, i).*qd(i));
                if i == 3
                    twistAxle = twist(:, i);
                end
            end
            
            % backward propagation of the articulated inertias and biases
            tau = nan(6, 1, class(q));
            tau(1:3) = 0;
            tau(4) = -obj.K(1)*(q(3)) - obj.C(1)*qd(4);  % vertical spring/damping
            tau(5) = -obj.K(2)*q(4) - obj.C(2)*qd(5);  % pitch torsional spring/damping
            tau(6) = -obj.K(3)*q(5) - obj.C(3)*qd(6);  % roll torsional spring/damping
            
            % inertia matrices
            zero = zeros(6, 6, class(q));
            II{1} = zero;
            II{2} = II{1};
            II{3} = obj.Maxle; % axle inertia
            II{4} = II{1};
            II{5} = II{1};
            II{6} = obj.M; % body inertia
            
            Mhat = cell(6, 1);             % init articulated inertia
            Mbar = cell(6, 1);             % init projected inertia
            bhat = nan(6, 6, class(q));    % init articulated bias
            Mscalar = nan(6, 1, class(q)); % init scalarazied non structural inertia (expression collection)
            
            for i = 6:-1:1
                if i == 6 % final body has no children
                    Mhat{i} = II{i};
                    Mscalar(i) = Xu(:, i).'*Mhat{i}*Xu(:, i) + 1e-3;
                    Mbar{i} = (  eye(6) - (Mhat{i}*Xu(:, i)*Xu(:, i).')./(Mscalar(i))  )*Mhat{i};
                    
                    bhat(:, i) = adStar(twist(:, i))*II{i}*twist(:, i);
                    g_axle_body = obj.Blocal{3}*obj.Blocal{4}*obj.Blocal{5};
                    rot_axle_body = g_axle_body(1:3,1:3);
                    rot_body_axle = rot_axle_body';
                    Faero_body = rot_body_axle*[-0.5*obj.cx*obj.rho*obj.Sf*twistAxle(1)^2,0,-0.5*(obj.cz1+obj.cz2)*obj.rho*obj.Sf*twistAxle(1)^2]';
                    Maero_body = rot_body_axle*[0,(obj.cz1*obj.a1-obj.cz2*obj.a2)*0.5*obj.rho*obj.Sf*twistAxle(1)^2,0]';
                    bhat(:, i) = bhat(:, i) - [Faero_body; Maero_body];
                    continue
                end
                
                j = i+1; % children body
                AdGstar = AdGinv{j}.';
                
                % current articulated inertia
                Mhat{i} = II{i} + AdGstar*Mbar{j}*AdGinv{j};
                Mscalar(i) = Xu(:, i).'*Mhat{i}*Xu(:, i) + 1e-3;
                
                % projected inertia (not required for the 1st body)
                if i ~= 1
                    Mbar{i} = (  eye(6) - (Mhat{i}*Xu(:, i)*Xu(:, i).')./(Mscalar(i))  )*Mhat{i};
                end
                
                % current body bias
                b_i = adStar(twist(:, i))*II{i}*twist(:, i);
                
                % current articulation bias
                V_ij = AdGinv{j}*twist(:, i); % current body twist in children frame
                bhat(:, i) = b_i + AdGstar*...
                    (...
                    bhat(:, j) - Mbar{j}*a{j}*V_ij + Mhat{j}*Xu(:, j).*(tau(j) - Xu(:, j).'*bhat(:, j))./Mscalar(j)...
                    ); % current articulated bias
                
                if i == 3 % axle body
                    bhat(:, i) = bhat(:, i) - wp;  % add pneumatic wrench
                end
            end
            
            % forward propagation of accelerations
            qdd = nan(6, 1, class(q));
            for i = 1:6
                if i == 1 % track joint transform
                    Vdot_parent = AdGinv{i}*baseAcc;
                    jac_acc = Jac_der_num.*alphadot.^2;
                    qdd(i) = -Xu(:, i).'*(Mhat{i}*(Vdot_parent + jac_acc) + bhat(:, i))./(Mscalar(i));
                    twistD(:, i) = Vdot_parent + Xu(:, i).*qdd(i) + jac_acc; % base velocity is always zero
                    continue
                end
                % current joint acceleration
                Vdot_parent = AdGinv{i}*twistD(:, i-1);
                V_p = a{i}*AdGinv{i}*twist(:, i-1);
                qdd(i) = (tau(i) - Xu(:, i).'*(Mhat{i}*(Vdot_parent - V_p) + bhat(:, i)))./(Mscalar(i));
                
                % current body twist derivative
                twistD(:, i) = Vdot_parent + Xu(:, i).*qdd(i) - V_p;
                
                % extract wrench on axles
                if i == 3
                    wrenchAxle = Mhat{i}*twistD(:, i) + bhat(:, i);
                end
            end
            
        end
        
        function [qdd, twistAxle, wrenchAxle] = ABAdyn_subGraph(obj, q, qd, wp, alpha_num, deltaAlpha) % tau is implicitly in the expression tau = fun(q, qd)
            %
            % [qdd, twistAxle, wrenchAxle] = ABAdyn(obj, Vbase, q, qd, wp)
            %
            % returns the virtual joint accelerations q_dotdotof the model
            % obtained with the Articulated Body Algorithm (ABA)
            %
            % states q = [alpha, ep, epsi, d, theta, phi]
            %        qdot = [alphadot, epdot, epsidot, ddot, thetadot, phidot]
            %
            % model input: wp: equivalent wrench applied to the axle frame
            %                 (third kinematic body)
            %
            % base twist is assumed null and the base acceleration is
            % assumed with an upward (along z) gravitational acceleration
            
            
            alpha = q(1);
            alphadot = qd(1);
            import casadi.*
            alpha = SX.sym('alpha');
            
            
            % update the track computational graph based on the knot span
            obj.track.update_subGraphBases(alpha, alpha_num, deltaAlpha);
            
            % compute some secondary functions
            % track body jacobian

            J = obj.track.fun_Vgs_s(alpha)*obj.track.nts;
            J(2) = 0; % simplify the computational graph by setting obvious zero values
            J(3) = 0;
            obj.T_jac = Function('T_jac', {alpha}, {J});
            
            % track body jacobian derivative
            Jder = jacobian(obj.track.fun_Vgs_s(alpha), alpha)*obj.track.nts + obj.track.fun_Vgs_s(alpha)*jacobian(obj.track.nts, alpha);
            Jder(2) = 0; % simplify the computational graph by setting obvious zero values
            Jder(3) = 0;
            obj.T_jac_der = Function('T_der', {alpha}, {Jder});
            
            alpha = alpha_num;
            
            % transform matrices update
            obj.FWKin_G(alpha, q); %%q(2:5)
            % add track jacobian as an equivalent unitary twist
            Xu = [full(obj.T_jac(alpha)), obj.X];
            
            baseAcc = [0;0;9.81;0;0;0]; % simulating an upward gravitational acceleration
            
            % twist forward propagation
            twist = nan(6, 6, class(q));  % twist init
            twistD = nan(6, 6, class(q)); % twist derivative initialization
            AdGinv = cell(1, 6); % adjoint matrices init
            a = cell(6, 1);                % init relative twist ad
            
            % first joint (track curvilinear joint) is a special case
            twist(:, 1) = Xu(:, 1).*alphadot; % twist of the track frame
            AdGinv{1} = adjointInv(obj.S);
            %            a{1} = ad(Xu(:, 1).*alphadot);
            
            % propagation loop of twits
            for i = 2:6
                AdGinv{i} = adjointInv(obj.Blocal{i-1});
                
                twist(:, i) = AdGinv{i}*twist(:, i-1) + Xu(:, i).*qd(i); % twist of i-th body (i = 1, is the base)
                a{i} = ad(Xu(:, i).*qd(i));
                if i == 3
                    twistAxle = twist(:, i);
                end
            end
            
            % backward propagation of the articulated inertias and biases
            tau = nan(6, 1, class(q));
            tau(1:3) = 0;
            tau(4) = -obj.K(1)*(q(3)) - obj.C(1)*qd(4);  % vertical spring/damping
            tau(5) = -obj.K(2)*q(4) - obj.C(2)*qd(5);  % pitch torsional spring/damping
            tau(6) = -obj.K(3)*q(5) - obj.C(3)*qd(6);  % roll torsional spring/damping
            
            % inertia matrices
            if isa(q, 'casadi.SX'); zero = casadi.SX(6,6); else; zero = zeros(6,6); end
            II{1} = zero;
            II{2} = II{1};
            II{3} = obj.Maxle; % axle inertia
            II{4} = II{1};
            II{5} = II{1};
            II{6} = obj.M; % body inertia
            
            Mhat = cell(6, 1);             % init articulated inertia
            Mbar = cell(6, 1);             % init projected inertia
            bhat = nan(6, 6, class(q));    % init articulated bias
            Mscalar = nan(6, 1, class(q)); % init scalarazied non structural inertia (expression collection)
            
            for i = 6:-1:1
                if i == 6 % final body has no children
                    Mhat{i} = II{i};
                    Mscalar(i) = Xu(:, i).'*Mhat{i}*Xu(:, i);
                    Mbar{i} = (  eye(6) - (Mhat{i}*Xu(:, i)*Xu(:, i).')./(Mscalar(i))  )*Mhat{i};
                    
                    bhat(:, i) = adStar(twist(:, i))*II{i}*twist(:, i);
                    g_axle_body = obj.Blocal{3}*obj.Blocal{4}*obj.Blocal{5};
                    rot_axle_body = g_axle_body(1:3,1:3);
                    rot_body_axle = rot_axle_body';
                    Faero_body = rot_body_axle*[-0.5*obj.cx*obj.rho*obj.Sf*twistAxle(1)^2,0,-0.5*(obj.cz1+obj.cz2)*obj.rho*obj.Sf*twistAxle(1)^2]';
                    Maero_body = rot_body_axle*[0,(obj.cz1*obj.a1-obj.cz2*obj.a2)*0.5*obj.rho*obj.Sf*twistAxle(1)^2,0]';
                    bhat(:, i) = bhat(:, i) - [Faero_body; Maero_body];
                    continue
                end
                
                j = i+1; % children body
                AdGstar = AdGinv{j}.';
                
                % current articulated inertia
                Mhat{i} = II{i} + AdGstar*Mbar{j}*AdGinv{j};
                Mscalar(i) = Xu(:, i).'*Mhat{i}*Xu(:, i);
                
                % projected inertia (not required for the 1st body)
                if i ~= 1
                    Mbar{i} = (  eye(6) - (Mhat{i}*Xu(:, i)*Xu(:, i).')./(Mscalar(i))  )*Mhat{i};
                end
                
                % current body bias
                b_i = adStar(twist(:, i))*II{i}*twist(:, i);
                
                % current articulation bias
                V_ij = AdGinv{j}*twist(:, i); % current body twist in children frame
                bhat(:, i) = b_i + AdGstar*...
                    (...
                    bhat(:, j) - Mbar{j}*a{j}*V_ij + Mhat{j}*Xu(:, j).*(tau(j) - Xu(:, j).'*bhat(:, j))./Mscalar(j)...
                    ); % current articulated bias
                
                if i == 3 % axle body
                    bhat(:, i) = bhat(:, i) - wp;  % add pneumatic wrench
                end
            end
            
            % forward propagation of accelerations
            qdd = nan(6, 1, class(q));
            for i = 1:6
                if i == 1 % track joint transform
                    Vdot_parent = AdGinv{i}*baseAcc;
                    jac_acc = obj.T_jac_der(alpha).*alphadot.^2;
                    qdd(i) = -Xu(:, i).'*(Mhat{i}*(Vdot_parent + jac_acc) + bhat(:, i))./(Mscalar(i));
                    twistD(:, i) = Vdot_parent + Xu(:, i).*qdd(i) + jac_acc; % base velocity is always zero
                    continue
                end
                % current joint acceleration
                Vdot_parent = AdGinv{i}*twistD(:, i-1);
                V_p = a{i}*AdGinv{i}*twist(:, i-1);
                qdd(i) = (tau(i) - Xu(:, i).'*(Mhat{i}*(Vdot_parent - V_p) + bhat(:, i)))./(Mscalar(i));
                
                % current body twist derivative
                twistD(:, i) = Vdot_parent + Xu(:, i).*qdd(i) - V_p;
                
                % extract wrench on axles
                if i == 3
                    wrenchAxle = Mhat{i}*twistD(:, i) + bhat(:, i);
                end
            end
            
        end
        
    end
    
    methods % visible methods
        
        function obj = track_vechicle_model(track, q0, options) % !!!***other inputs will be added when more functionality is implemented (inertial properties, springs, dampers, tires etc***
            % 
            % track is an object or struct with fields: 
            % x : x(s) position of the mid-line track curve
            % t : t(s) tangent vector
            % n : n(s) transversal vector
            % m : m(s) normal vector (normal to the road)
            % 
            % Vgs : Vgs(s) twist ( motion of {S} w.r.t. {G} (ground) as
            %       function of s, with components written in {S}
            %
            % q0 is the initial configuration of the vehicle w.r.t. the
            % track: [e_p0; e_psi0; d0; theta0; phi0]
            %
            arguments
               
                track
                q0
                options.K      = [];
                options.C      = [];
                options.a1     = [];
                options.a2     = [];
                options.cz1    = [];
                options.cz2    = [];
                options.cx     = [];
                options.rho    = [];
                options.Sf     = [];
                options.q1     = [];
                options.q2     = [];
                options.qm     = [];
                options.h      = [];
                options.GM     = [];
                options.k_phi1 = [];
                options.k_phi2 = [];
                options.k_phi = [];
                options.Mbody  = [];
                options.Maxle  = [];
                
            end
            obj.K      = options.K;
            obj.C      = options.C;
            obj.a1     = options.a1;
            obj.a2     = options.a2;
            obj.cz1    = options.cz1;
            obj.cz2    = options.cz2;
            obj.cx     = options.cx;
            obj.rho    = options.rho;
            obj.Sf     = options.Sf;
            obj.q1     = options.q1;
            obj.q2     = options.q2;
            obj.qm     = options.qm;
            obj.h      = options.h;
            obj.GM     = options.GM;
            obj.k_phi1 = options.k_phi1;
            obj.k_phi2 = options.k_phi2;
            obj.k_phi  = options.k_phi;
            obj.M      = options.Mbody;
            obj.Maxle  = options.Maxle;
            
            % extract and save inputs
            obj.track = track;
            obj.q0 = q0;
            
            % define the unitary twists once
            obj.X =          [[0;1;0;0;0;0],... % e_p
                              [0;0;0;0;0;1],... % e_psi
                              [0;0;1;0;0;0],... % d
                              [0;0;0;0;1;0],... % theta
                              [0;0;0;1;0;0]];   % phi
                          
            % define the offset transform matrices once
            obj.g0{1} = Tty(q0(1));
            obj.g0{2} = TrotZ(q0(2));
            obj.g0{3} = Ttz(q0(3));
            obj.g0{4} = TrotY(q0(4));
            obj.g0{5} = TrotX(q0(5));
            
            % compute some secondary functions
            import casadi.*
            x = MX.sym('x');
            % track body jacobian
            J = track.fun_Vgs_s(x)*track.fun_normts(x);
            J(2) = 0; % simplify the computational graph by setting obvious zero values
            J(3) = 0;
            obj.T_jac = Function('T', {x}, {J});
            % track body jacobian derivative
            Jder = jacobian(track.fun_Vgs_s(x), x)*track.fun_normts(x) + track.fun_Vgs_s(x)*jacobian(track.fun_normts(x), x);
            Jder(2) = 0; % simplify the computational graph by setting obvious zero values
            Jder(3) = 0;
            obj.T_jac_der = Function('T_der', {x}, {Jder});
        end
        
        function J = ground_to_axle_jacobian(obj, q)
            %
            % returns the body jacobian matrix from ground to axle (axle
            % frame)
            %
            % q = [alpha, ep, ef, d, theta, phi]
            
            [~, T] = obj.FWKin(q(2:end));
            A_sa = T{1}*T{2};
            
            Js = obj.baseTwist(q(1), 0);
            J = [adjoint(rigidInverse(A_sa))*Js, obj.track_to_axle_jacobian(q(2:end))];
        end
        
        function [B, Blocal] = FWKin(obj, q)
            %
            % FWKin(this, q)
            % forward kinematics computed with the local POE formulation
            % returns the transform matrix from track frame {S} to body {B}
            % 
            if nargout > 1
                [B, Blocal] = FWkin_localPOE(obj.g0, obj.X, q);
                return
            end
            [obj.B, obj.Blocal] = FWkin_localPOE(obj.g0, obj.X, q);

            
            % since the unitary twists can coincide with the global unitary
            % twists, the kinematics can be formulated with the global POE
            % at a specific initial configuration ({B} coincides with {S}),
            % but we lose information about the local transforms
            
            %             obj.B = FWKin(eye(4), {obj.X(:,1), obj.q0(1) + q(1)},...
            %                 {obj.X(:,2), obj.q0(2) + q(2)},...
            %                 {obj.X(:,3), obj.q0(3) + q(3)},...
            %                 {obj.X(:,4), obj.q0(4) + q(4)},...
            %                 {obj.X(:,5), obj.q0(5) + q(5)});

        end
        
        function FWKin_G(obj, alpha, q)
            %
            % forward kinematics from ground {G} to body {B}
            %
            
            obj.FWKin(q);                            % kinem from {S} to body {B}
            obj.S = full(obj.track.fun_g_gs(alpha)); % kinem from {G} to body {S}
            obj.Ggb = obj.S*obj.B;                   %
            
        end
        
        function [xdot, w2, V_g2] = state_space_eq(obj, alpha, q, V_b, wp)
            % 
            % INPUTS:
            %   alpha : track variable
            %   q: kinematic configuration variables [ep, epsi, d, theta, phi]
            %   Vb2: body {B} w.r.t. ground {G} twist written in {2} (axle)
            %   wp: penumatic wrench given in {2}
            %
            % OUTPUTS:
            %   xdot: state derivative vector : [alphadot; qdot; Vb2dot]
            %         12x1
            %   w2: wrench applied to the axle in {2}
            
            % define a zero 3x3 matrix with structural zeros if casadi is used
            if isa(q, 'casadi.SX'); import casadi.*; zero = SX(3,3); else; zero = zeros(3,3); end
            
            % kinematic computation to update the transforms
            obj.FWKin_G(alpha, q);
            
            % base twist computation
            [~, Jac_base] = obj.baseTwist(alpha, 1);
            
            % ground to body jacobian computation in {B}
            G = obj.Blocal{3}*obj.Blocal{4}*obj.Blocal{5};
            G_sb = obj.B;
            J = [adjoint(rigidInverse(G_sb))*full(Jac_base), obj.track_to_body_jacobian(q)];
            
            % positional states dynamic
            qd = J\V_b; % = [alpha_dot, e_p_dot, e_psi_dot, d_dot, theta_dot, phi_dot]^T
            
            % gravity wrench
            wg_sB = [0;0;-9.81;0;0;0].*obj.M(1,1); % gravity wrench in {S} applied on Ob (does not retain the Lie parametrization convention)
            Rgb = obj.Ggb(1:3, 1:3); 
            wg_b = [Rgb.', zero; zero, zero]*wg_sB; % gravity wrench in {B} applied on Ob
            
            % dynamic equilibrium of the body
            Vdot_b = obj.bodyDyn(q(3:5), V_b, wp, wg_b, qd(4:6));
            
            % state-space dynamics
            xdot = [qd; Vdot_b];
            
            % computation of structural wrench in {2}
            wb = -wg_b + obj.M*Vdot_b + adStar(V_b)*obj.M*V_b;
            w2 = adjointStar(G)*wb;
            
            % ground to Axle jacobian
            G_s2 = obj.Blocal{1}*obj.Blocal{2};
            J_g2 = [adjoint(rigidInverse(G_s2))*full(Jac_base), obj.track_to_axle_jacobian(q)];
            V_g2 = J_g2*[qd(1); qd(2); qd(3)];
        end
        
        function computeCasadi_state_space_function(obj, options)
           %
           % 
           %
           arguments 
               obj
               options.compile = false;
               options.integration_var = 'time';
           end
           import casadi.*
           % init variables
           qsym = MX.sym('q', 6, 1);
           qdotsym = MX.sym('qdot', 6, 1);           
           wsym = MX.sym('w', 6, 1);
           
           % dynamics
           [qdd, V_g2, wA] = obj.ABAdyn(qsym, qdotsym, wsym);
           
           % state space field
           xdot = vertcat(qdotsym, qdd);
           
           if ~strcmpi(options.integration_var, 'time')
               xdot = xdot./xdot(1); % normalizing w.r.t. alpha if the integration variable is not the time
           end
           
           % function creation
           obj.xdot_sym = Function('sdot', {[qsym; qdotsym], wsym}, {xdot});
           obj.axle_wrench_sym = Function('wAxle', {[qsym; qdotsym], wsym}, {wA});
           obj.axle_twist = Function('Vaxle', {[qsym; qdotsym]}, {V_g2});
           
           if options.compile 
               options = struct('mex', true, 'main', false);
               obj.xdot_sym.generate('xdot_car.c', options);
               obj.axle_wrench_sym.generate('wrench_car.c', options);
               obj.axle_twist.generate('twist_car.c', options);

               mex xdot_car.c -largeArrayDims
               mex wrench_car.c -largeArrayDims
               mex twist_car.c -largeArrayDims
           end 
        end
        
        function computeCasadi_state_space_function_alpha_num(obj, options)
            %
            %
            %
            arguments
                obj
                options.compile = false;
                options.integration_var = 'time';
            end
            import casadi.*
            % init variables
            g_gs_sym = SX.sym('g_gs', 4, 4);
            Jac_sym = SX.sym('Jac', 6, 1);
            Jac_der_sym = SX.sym('Jacdot', 6, 1);
            qsym = SX.sym('q', 5, 1);
            qdotsym = SX.sym('qdot', 6, 1);
            wsym = SX.sym('w', 6, 1);
            
            % dynamics
            [qdd, V_g2, wA] = obj.ABAdyn_alpha_num(g_gs_sym, Jac_sym, Jac_der_sym, qsym, qdotsym, wsym);
            
            % state space field
            xdot = vertcat(qdotsym, qdd);
            
            if ~strcmpi(options.integration_var, 'time')
                xdot = xdot./xdot(1); % normalizing w.r.t. alpha if the integration variable is not the time
            end
            
            % function creation
            obj.xdot_sym = Function('sdot', {g_gs_sym, Jac_sym, Jac_der_sym, [qsym; qdotsym], wsym}, {xdot});
            obj.axle_wrench_sym = Function('wAxle', {g_gs_sym, Jac_sym, Jac_der_sym, [qsym; qdotsym], wsym}, {wA});
            obj.axle_twist = Function('Vaxle', {g_gs_sym, Jac_sym, Jac_der_sym, [qsym; qdotsym]}, {V_g2});
            
            if options.compile
                options = struct('mex', true, 'main', false);
                obj.xdot_sym.generate('xdot_car.c', options);
                obj.axle_wrench_sym.generate('wrench_car.c', options);
                obj.axle_twist.generate('twist_car.c', options);
                
                mex xdot_car.c -largeArrayDims
                mex wrench_car.c -largeArrayDims
                mex twist_car.c -largeArrayDims
            end
        end
        
        function out = computeRK4step(obj, q, qd, wp, dt)
           %
           %
           %
           % 
           if isempty(obj.xdot_sym)
               obj.computeCasadi_state_space_function;
           end 
           
           % Penumatic constitutive equation
           % controls should replace wp in the inputs
           % w_p = ...
           
           % dynamic equations
           %[snext, obj.BS_solver] = Euler_BS([q; qd], wp, @obj.xdot_sym, dt, 'Solver', obj.BS_solver); % RK4_step([q; qd], wp, @xdot_car, dt, 'Nsteps', 2);
           snext = RK4_step([q; qd], wp, @obj.xdot_sym, dt);
           snext = full(snext);
           
           if isnan(snext(1))
               error('NaN detected: reached the end of the track, alpha > 1 or unstable integration')
           end
           
           % extract and update states
           obj.q = snext(1:6); % configuration
           obj.qdot = snext(7:12); % configuration derivatives
           
           % axle wrench
           obj.w2 = full(obj.axle_wrench_sym([q; qd], wp));
           
           % axle twist
           obj.V = full(obj.axle_twist([q; qd]));
           
           obj.FWKin_G(obj.q(1), obj.q(2:6)); % update transform matrices (this is just required for the real-time plotting)
           
           % output if prompted
           if nargout > 0
              out = snext; 
           end
        end
        
        function out = compute_Direct_Collocation_step(obj, q, qd, wp, dt, options)
           %
           %
           %
           % 
           arguments
                obj
                q
                qd
                wp
                dt
                options.degree = 2;           
           end
           
           d = options.degree;
           if isempty(obj.xdot_sym)
               obj.computeCasadi_state_space_function;
           end 
           
           % Penumatic constitutive equation
           % controls should replace wp in the inputs
           % w_p = ...

           % dynamic equations
           snext = Direct_collocation_num([q; qd], wp, @obj.xdot_sym, dt, 'degree', d);
           snext = full(snext);
           
           if isnan(snext(1))
               error('NaN detected: reached the end of the track, alpha > 1 or unstable integration')
           end
           
           % extract and update states
           obj.q = snext(1:6); % configuration
           obj.qdot = snext(7:12); % configuration derivatives
           
           % axle wrench
           obj.w2 = full(obj.axle_wrench_sym([q; qd], wp));
           %obj.w2 = full(obj.axle_wrench_sym([q; qd], wp));
           % axle twist
           obj.V = full(obj.axle_twist([q; qd]));
           %obj.V = full(obj.axle_twist([q; qd]));
           
           obj.FWKin_G(obj.q(1), obj.q(2:6)); % update transform matrices (this is just required for the real-time plotting)
           
           % output if prompted
           if nargout > 0
              out = snext; 
           end
        end
        
        function ax = plot(obj, options)
            %
            % plots the vehicle on the track
            %
            
            arguments
                obj
                options.plot_track = true;
                options.alpha_end = 1;
                options.parent = [];
                options.labels = true;
                
            end
            
            % import the stl placeholder model
            [Faces, Vertices] = stlread('F1_2016_car.stl');
            %adjust vertices
            %Vertices = Vertices - [+2.5, +1.7/2, 0]; %1.25/2 FOR EUG MESH
            Vertices = Vertices*0.015;
            Vertices = (rotZ(pi/2)*Vertices')';
            
            if isempty(options.parent)
                figure; hold on; axis equal; box on; grid on;
                parent = gca;
                
            else
                parent = options.parent;
            end
            % initialize transform for spatial frame
            obj.transformS = hgtransform(parent);
            % initialize transform for body frame
            obj.transformB = hgtransform(obj.transformS);
            obj.transformB.Matrix = obj.B; % from {S} to {B}
            obj.transformS.Matrix = obj.S; % from {G} to {S}
            
            % car patch plot
            patch('Faces', Faces, 'Vertices', Vertices, 'edgecolor', 'none', 'faceColor', 'r', 'parent', obj.transformB);
            
            if options.plot_track
                % track plot
                obj.track.full_plot(linspace(0,options.alpha_end,1500), 'parent', parent);
            end
            
            % frames plot
            if options.labels
                lab = {'s', 'b'};
            else
                lab = {[], []};
            end
            plotFrame(eye(4), 'parent', obj.transformS, 'label', lab{1})
            plotFrame(eye(4), 'parent', obj.transformB, 'label', lab{1})
            
            % axis settings
            view(parent, 45,25)
            camlight(parent, 'headlight', 'infinite')
            lighting(parent, 'gouraud')

            obj.axes = parent;
            camproj(obj.axes, 'perspective')
            camva(obj.axes, 35)
            
            if nargout > 0 % return the plot axes
               ax = obj.axes; 
            end
        end
        
        function updateplot(obj, options)
            %
            % refreshes the plot with updated states
            % OPTIONS:
            %  - camera: 1 (top view); 0 (3d view)
            
            arguments
               obj
               options.camera = 0;
               options.boxSize = [20, 20, 20];
            end
            obj.transformB.Matrix = obj.B;
            obj.transformS.Matrix = obj.S;
            G = obj.Ggb;
            pos = G(1:3, 4);
            % center the plot around the computed track point
            if options.camera == 1
                t = G(1:3, 1);
                w = G(1:3, 3);
                tr = obj.S(1:3, 1);
                wtr = obj.S(1:3, 3);
                campos(obj.axes, G(1:3, 4) - t*50 + w*25)
                camtarget(obj.axes, G(1:3, 4) + t*15 - w*2)
            else
                set(obj.axes, 'xlim', [pos(1)-options.boxSize(1),  pos(1)+options.boxSize(1)])
                set(obj.axes, 'ylim', [pos(2)-options.boxSize(2),  pos(2)+options.boxSize(2)])
                set(obj.axes, 'zlim', [pos(3)-options.boxSize(3),  pos(3)+options.boxSize(3)])
                
            end

        end
        
        function update_state_space_functions(obj, alpha_num, deltaAlpha, options)
           %
           % 
           %
           arguments 
               obj
               alpha_num
               deltaAlpha
               options.compile = false;
               options.integration_var = 'alpha';
           end
           import casadi.*
           % init variables
           qsym = SX.sym('q', 5, 1);
           qdotsym = SX.sym('qdot', 6, 1);           
           wsym = SX.sym('w', 6, 1);
           
           % dynamics
           [qdd, V_g2, wA] = obj.ABAdyn_subGraph(qsym, qdotsym, wsym, alpha_num, deltaAlpha);
           
           % state space field
           xdot = vertcat(qdotsym, qdd);
           
           if ~strcmpi(options.integration_var, 'time')
               xdot = xdot./xdot(1); % normalizing w.r.t. alpha if the integration variable is not the time
           end
           
           % function creation
           obj.xdot_sym = Function('sdot', {[qsym; qdotsym], wsym}, {xdot});
           obj.axle_wrench_sym = Function('wAxle', {[qsym; qdotsym], wsym}, {wA});
           obj.axle_twist = Function('Vaxle', {[qsym; qdotsym]}, {V_g2});
           
           if options.compile 
               options = struct('mex', true, 'main', false);
               obj.xdot_sym.generate('xdot_car.c', options);
               obj.axle_wrench_sym.generate('wrench_car.c', options);
               obj.axle_twist.generate('twist_car.c', options);

               mex xdot_car.c -largeArrayDims
               mex wrench_car.c -largeArrayDims
               mex twist_car.c -largeArrayDims
           end 
        end
        
    end
    

end