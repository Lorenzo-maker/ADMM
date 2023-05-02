classdef nurbs_ribbon_casadi < handle % lead to change class property with the execution of a function 
        
    properties
    x               % x coordinates of the track
    y               % y coordinates of the track
    z               % z coordinates of the track 
    wl              % left width
    wr              % right width
    tw              % twist angle 
    n               % number of control points -1
    p               % degree of spline (oder = p + 1)
    m               % number of knots point -1
    in              % index of starting point
    out             % index of ending point
    ub              % vector of equations points
    u               % vector of knot points
    data            % structur with some other useful data 
    todo            % 'interp' or 'fit'
    q               % number of track point -1. Useful for knots point in 'fit' case
    N               % basis function
    Np              % last column of basis function 
    Nd              % basis of derivative
    Ndp             % last column of derivatives basis
    Ndd             % basis of second derivative
    Nddp            % last column of second derivatives basis
    ts              % tangent vector
    nts             % tangent vector norm
    nv              % normal vector (from derivative) 
    vh              % normal vector (from rotation) 
    k               % theorical curvature (always positive)
    kv              % theorical vector of curvature
    kp              % curvature in vh direction  
    kns             % curvature in ns direction
    bnv             % binormal vector (from second derivative)
    wh              % binormal vector (from rotate tangent)
    ubm             % method to build ubar
    knm             % method to build knots
    A               % matrix of interp problem
    Px              % x-coordinates of control points
    Py              % y-coordinates of control points
    Pz              % z-coordinates of control points 
    Pwl             % control points for left width nurbs
    Pwr             % control points for right width nurbs
    Ptw             % control points for twist nurbs
    q_error         % vector of quadratic error
    tot_error       % total quadratc error
    NN              % matrix for fit problem
    limited         % boolean for nurbs limitation to alfa_lim
    alfa_lim        % [alfa_init, alfa_end] for nurbs limitation
    Points          % matrix with x-y-z of external, middle and internal line
    fun_basis       % casadi function to compute symbolic and numerical basis
    fun_basis_Np    % casadi function to compute last column of basis matrix N 
    f_example       % casadi function to test basis, given p, us and knots point
    fun_pos         % casadi function to obtain numerical or symbolic position from us
    fun_ts          % casadi function to evaluate tangent vector of S-frame write in G-frame 
    fun_normts      % casadi function to evaluate the norm of tangent vector, symbolic or numerical
    fun_nv          % casadi function to evaluate components of normal vector (from derivative), symbolic or numerical
    fun_vh          % casadi function to evaluate normal vector of H-frame write in G-frame
    fun_bnv         % casadi function to evaluate components of binormal vector (from second derivative), symbolic or numerical
    fun_wh          % casadi function to evaluate components of binormal vector (from rotate tangent), symbolic or numerical
    fun_k           % casadi function to evaluate theorical curvature, symbolic or numerical 
    fun_kv          % casadi function to evaluate theorical curvature vector of binormal vector, symbolic or numerical
    fun_kp          % casadi function to evaluate curvature from vh(ex pt), symbolic or numerical
	fun_Rgh         % casadi function to evaluate rotation matrix from H-frame to S-frame
    fun_Rgs         % casadi function to evaluate rotation matrix from S-frame to G-frame
    fun_ns          % casadi function to evaluate normal vector of S-frame write in G-frame
    fun_ms          % casadi function to evaluate binormal vector of S-frame write in G-frame
    fun_g_gs        % casadi function to evaluate vector transformation from S-frame to G-frame
    fun_kns         % casadi function to evaluate curvature in ns direction
    fun_Vgs_s       % casadi function to evaluate twist geometric vector of S respect to G write in S frame
    fun_mudot
    end
        
    methods 
        %% method to inizialize the class
        function obj = nurbs_ribbon_casadi(track_point, in, out, p, n, todo, ubar_method, knots_method, alfa_limited, alfa_lim) 
            %obj = threeD_nurbs_problem_casadi(track_point, in, out, p, n, todo, ubar_method, knots_method) 
            %
            % matlab class necessary to build a NURBS model of the track
            %
            % track_point: track points that have to be described by a NURBS
            %
            % in: enter the starting point of the track
            % 
            % out: enter the ending point of the track (enter 'end' to get to the end)
            %
            % p: the spline degree
            %
            % n: number of spline control points
            %
            % todo: the method to approximate the spline 
            %        - 'interp' : to interpolate the track with the spline
            %        - 'fit' : to fit the track with the spline 
            %
            % ubar_method: the method to calculate the point for imposing equalities (Nurbs book pag.364-365)
            %               - chord 
            %               - centripetal
            %
            % knots_metod: the method to calculate the knot points (Nurbs book pag.365/pag.412)
            %               - equally (available only for interpolating problem)
            %               - averaging (available only for interpolating problem)
            %               - deboor (available only for fitting problem) 
            obj.limited = alfa_limited;
            obj.alfa_lim = alfa_lim;
            obj.ubm = ubar_method;
            obj.knm = knots_method;
            obj.in = in;
            dim0 = length(track_point(:,1));
            if strcmp(out, 'end') 
                obj.out = dim0;
            else
                obj.out = out;
            end
            obj.x  = track_point(obj.in:obj.out,1);
            obj.y  = track_point(obj.in:obj.out,2);
            obj.z  = track_point(obj.in:obj.out,3);
            obj.wl = track_point(obj.in:obj.out,4);
            obj.wr = track_point(obj.in:obj.out,5);
            obj.tw = track_point(obj.in:obj.out,6);
            obj.p = p;
            obj.q = length(obj.x)-1;
            obj.todo = todo;
            if strcmp(obj.todo, 'interp')
                obj.n = length(obj.x)-1;
            elseif strcmp(obj.todo, 'fit')
                if n < length(obj.x)-1
                    obj.n = n;
                else
                    error('fit is not possible, number of ctrl points exceed number of track point')
                end
            end
            obj.m = obj.p + obj.n + 1;
            obj.data.sqr = zeros(length(obj.x),1);
            for j = 2:length(obj.x)
                obj.data.sqr(j) = sqrt((obj.x(j)-obj.x(j-1)).^2 + (obj.y(j)-obj.y(j-1)).^2 + (obj.z(j)-obj.z(j-1)).^2);
            end
        end
        %% method to build eval points (pag.364-365 Nurbs book)
        function  chord(obj) 
            % method to build eval points (pag.364-365 Nurbs book)
            obj.data.d = sum(obj.data.sqr); 
            for j = 1:length(obj.x)
                obj.ub(j) = sum(obj.data.sqr(1:j));
            end
            obj.ub = obj.ub/obj.data.d;
        end
        function  centripetal(obj)
            % method to build eval points (pag.364-365 Nurbs book)
            obj.data.d = sum(sqrt(obj.data.sqr)); 
            for j = 1:length(obj.x)
                obj.ub(j) = sum(sqrt(obj.data.sqr(1:j)));
            end
            obj.ub = obj.ub/obj.data.d;
        end
        %% method to build knot points (pag.365/pag.412 Nurbs book)
        function equally(obj)
            % method to build knot points see pag.365 of NURBS book
            if strcmp(obj.todo, 'interp')
                obj.u = zeros(obj.p + 1, 1);
                for j = 1:(obj.n-obj.p) 
                    obj.u(j+obj.p+1) = j/(obj.n - obj.p + 1);
                end
                obj.u(obj.m - obj.p +1:obj.m + 1) = ones(obj.p + 1,1);
            else
                error('error: this method is not available for curve fitting')
            end
        end
        function averaging(obj)
            % method to build knot points see pag.365 of NURBS book
            if strcmp(obj.todo, 'interp')
                obj.u = zeros(obj.p + 1, 1);
                for j = 1:(obj.n-obj.p) 
                    obj.u(j+obj.p+1) = sum(obj.ub(j+1:j+obj.p))/obj.p;
                end
                obj.u(obj.m - obj.p +1:obj.m + 1) = ones(obj.p + 1,1);
            else 
                error('error: this method is not available for curve fitting')
            end
        end
        function deboor(obj)
            % method to build knot points see pag.412 of NURBS book
            if strcmp(obj.todo, 'fit')
                obj.u = zeros(obj.p + 1, 1);
                d = (obj.q + 1)/(obj.n - obj.p + 1);
                for j = 1:(obj.n-obj.p)
                    i = floor(j*d);
                    alfa = j*d - i;
                    obj.u(j+obj.p+1) = (1-alfa)*obj.ub(i)+alfa*obj.ub(i+1);
                end
                obj.u(obj.m - obj.p +1:obj.m + 1) = ones(obj.p + 1,1);
            else 
                error('error: this method is not available for curve interp')
            end
        end
        %% Basis function
        function  basis(obj,us) 
           arguments 
               obj
               us
           end
           
            % method to build basis function for NURBS construction
            %fornire in input SX.sym('us')
            %build basis function with casadi structure, if Barcellona is
            %our instance we run Barcellona.basis(SX.sym('us')), then we
            %can use Barcellona.fun_basis/Barcellona.fun_basis_Np(us) with
            %us numerical or symbolic to obtain numerical or symbolic basis
            import casadi.*
            if obj.limited
                disp('Optimization of nurbs computational graph...')
                index_in = find(obj.u <= obj.alfa_lim(1));
                index_in = index_in(end);
                index_end = find(obj.u >= obj.alfa_lim(2));
                index_end = index_end(1);  
            end
            
            obj.N = SX(obj.n+2,obj.p+1);
            obj.Np = [];
            for i = 1:(obj.n+1)
                if obj.limited
                    if obj.u(i) > obj.u(index_end) || obj.u(i+1) < obj.u(index_in)
                        continue
                    end
                end
                    obj.N(i,1) = if_else(((obj.u(i)) <= us & us < (obj.u(i+1))),1,DM(1,1));
                    obj.N(i,1) = if_else((obj.u(i)) < 1.0 & us ==1.0 & (obj.u(i+1)) == 1,1,obj.N(i,1));																				   
                    obj.N(i,2:obj.p+1) = DM(1,1);
                
            end
            obj.N(obj.n+2,:) = DM(1,1);
            for j = 1:obj.p    
                for i = 1:(obj.n+1)
                    num1 = if_else(DM(obj.u(i+j) - obj.u(i)) == 0,0,us - obj.u(i));
                    den1 = if_else(DM(obj.u(i+j) - obj.u(i)) == 0,1, obj.u(i+j) - obj.u(i));
                    num2 = if_else(DM(obj.u(i+j+1) - obj.u(i+1)) == 0,0,obj.u(i+j+1) - us);
                    den2 = if_else(DM(obj.u(i+j+1) - obj.u(i+1)) == 0,1,obj.u(i+j+1) - obj.u(i+1));
                    obj.N(i,j+1) = obj.N(i,j)*(num1/den1)  + obj.N(i+1,j)*(num2/den2);
                end
                    
            end
                obj.Np = obj.N(1:obj.n+1,obj.p+1);
                obj.fun_basis = Function('fun_basis',{us},{obj.N},{'us'},{'N'});
                obj.fun_basis_Np = Function('fun_basis_Np',{us},{obj.Np},{'us'},{'Np'});
        end
        
        function Np = basis_numerical(obj,us)
            % numerical basis function given a value of us
            obj.N = zeros(obj.n+2,obj.p+1);
            for i = 1:(obj.n+1)
                if (obj.u(i) <= us) && (us < obj.u(i+1))
                    obj.N(i,1) = 1;
                else
                    obj.N(i,1) = 0;
                end
            end
            for j = 1:obj.p    
                for i = 1:(obj.n+1)    
                    if obj.u(i+j) - obj.u(i) == 0
                        num1 = 0;
                        den1 = 1;
                    else
                        num1 = us - obj.u(i);
                        den1 = obj.u(i+j) - obj.u(i);
                    end
                    if obj.u(i+j+1) - obj.u(i+1) == 0
                        num2 = 0;
                        den2 = 1;
                    else
                        num2 = obj.u(i+j+1) - us;
                        den2 = obj.u(i+j+1) - obj.u(i+1);
                    end
                    obj.N(i,j+1) = obj.N(i,j)*(num1/den1)  + obj.N(i+1,j)*(num2/den2);
                    Np = obj.N(1:obj.n+1,obj.p+1);
                    obj.Np = Np;
                end
            end          
        end
          
        function basis_example(obj,us,pt,knots)
            % function to evaluate numerical or symbolic basis given order p and knots vector
            %input us = SX.sym('us')
            import casadi.*
            mt = (length(knots)-1);
            nt = mt - pt - 1;
            Nt = SX.sym('Nt',nt+2,pt+1);
            for i = 1:(nt+1)
                Nt(i,1) = if_else((DM(knots(i)) <= us & us < DM(knots(i+1))),1,0);
                Nt(i,2:pt+1) = DM(0); 
            end 
            Nt(nt+2,:) = DM(0);
            for j = 1:pt
                for i = 1:(nt+1)
                    num1 = if_else(DM(knots(i+j) - knots(i)) == 0,0,us - knots(i));
                    den1 = if_else(DM(knots(i+j) - knots(i)) == 0,1, knots(i+j) - knots(i));
                    num2 = if_else(DM(knots(i+j+1) - knots(i+1)) == 0,0, knots(i+j+1) - us);
                    den2 = if_else(DM(knots(i+j+1) - knots(i+1)) == 0,1,knots(i+j+1) - knots(i+1));
                    Nt(i,j+1) = Nt(i,j)*(num1/den1)  + Nt(i+1,j)*(num2/den2);
                end
            end
            Ntp = Nt(:,pt+1);
            obj.f_example = Function('f_example',{us},{Nt,Ntp},{'us'},{'Nt','Ntp'});
                         
        end  
        %% problem solution
        function solve(obj)
            % obtain the solution of the numerical problem of interp or fit, in order to obtain the coordinates of Px, Py, Pz. No casadi structure is needed here, the problem is always numeric.
            if strcmp(obj.todo, 'interp')
                if strcmp(obj.ubm, 'chord')
                    obj.chord;
                elseif strcmp(obj.ubm, 'centripetal')    
                    obj.centripetal;
                end
                if strcmp(obj.knm, 'equally')
                    obj.equally
                elseif strcmp(obj.knm, 'averaging')
                    obj.averaging
                elseif strcmp(obj.knm, 'deboor') 
                    error('The deboor method for knots is not available for interpolation problem')
                end
                obj.A = zeros(obj.q+1, obj.n+1);
                for i = 1:obj.q+1
                    obj.A(i,:) = obj.basis_numerical(min(0.999999,obj.ub(i)))';
                end
                obj.Px  = linsolve(obj.A,obj.x);
                obj.Py  = linsolve(obj.A,obj.y);
                obj.Pz  = linsolve(obj.A,obj.z);
                obj.Pwl = linsolve(obj.A,obj.wl);
                obj.Pwr = linsolve(obj.A,obj.wr);
                obj.Ptw = linsolve(obj.A,obj.tw);
            elseif strcmp(obj.todo, 'fit')
                 if strcmp(obj.ubm, 'chord')
                    obj.chord;
                 elseif strcmp(obj.ubm, 'centripetal')    
                    obj.centripetal;
                 end
                 if strcmp(obj.knm, 'equally')
                    error('The equally method for knots is not available for fitting problem')
                 elseif strcmp(obj.knm, 'averaging')
                    error('The averaging method for knots is not available for interpolation problem')
                 elseif strcmp(obj.knm, 'deboor') 
                    obj.deboor
                 end
                 obj.NN = zeros(obj.q-1, obj.n-1);
                 Rx  = zeros(obj.q-1,1);
                 Ry  = zeros(obj.q-1,1);
                 Rz  = zeros(obj.q-1,1);
                 Rwl = zeros(obj.q-1,1);
                 Rwr = zeros(obj.q-1,1);
                 Rtw = zeros(obj.q-1,1);
                 for i = 1:obj.q-1
                    %Nb =  obj.basis_numerical(min(0.999999,obj.ub(i)))';
                    Nb =  obj.basis_numerical(obj.ub(i+1))';
                    obj.NN(i,:) = Nb(2:obj.n);
                    Rx(i)  = obj.x(i+1) - Nb(1)*obj.x(1) - Nb(end)*obj.x(end);
                    Ry(i)  = obj.y(i+1) - Nb(1)*obj.y(1) - Nb(end)*obj.y(end);
                    Rz(i)  = obj.z(i+1) - Nb(1)*obj.z(1) - Nb(end)*obj.z(end);
                    Rwl(i) = obj.wl(i+1) - Nb(1)*obj.wl(1) - Nb(end)*obj.wl(end);
                    Rwr(i) = obj.wr(i+1) - Nb(1)*obj.wr(1) - Nb(end)*obj.wr(end);
                    Rtw(i) = obj.tw(i+1) - Nb(1)*obj.tw(1) - Nb(end)*obj.tw(end);
                 end
                 Rx  = obj.NN'*Rx;
                 Ry  = obj.NN'*Ry;
                 Rz  = obj.NN'*Rz;
                 Rwl = obj.NN'*Rwl;
                 Rwr = obj.NN'*Rwr;
                 Rtw = obj.NN'*Rtw;
                 
                 obj.Px = [obj.x(1);linsolve(obj.NN'*obj.NN,Rx);obj.x(end)];
                 obj.Py = [obj.y(1);linsolve(obj.NN'*obj.NN,Ry);obj.y(end)];
                 obj.Pz = [obj.z(1);linsolve(obj.NN'*obj.NN,Rz);obj.z(end)];
                 obj.Pwl = [obj.wl(1);linsolve(obj.NN'*obj.NN,Rwl);obj.wl(end)];
                 obj.Pwr = [obj.wr(1);linsolve(obj.NN'*obj.NN,Rwr);obj.wr(end)];
                 obj.Ptw = [obj.tw(1);linsolve(obj.NN'*obj.NN,Rtw);obj.tw(end)];
            end
                
        end
        %% evaluate nurbs    
        function  eval(obj,us)
            % build the casadi function fun_pos(us) necessary to evaluate the spline
            % insert us = SX.sym('us')
            %if Barcellona is our instance we run
            %Barcellona.eval_symbolic(SX.sym('us')) and then we have
            %Barcellona.fun_pos(us) with us numerical or symbolic to obtain
            %symbolic or numerical position. us could be a vector. 
            import casadi.*
            if isempty(obj.Px) 
                error('First run obj.solve')
            end
            if isempty(obj.fun_basis_Np)
                obj.basis(us);
            end
            x_p = dot(obj.fun_basis_Np(us),DM(obj.Px));
            y_p = dot(obj.fun_basis_Np(us),DM(obj.Py));
            z_p = dot(obj.fun_basis_Np(us),DM(obj.Pz));
            wl_p = dot(obj.fun_basis_Np(us),DM(obj.Pwl));
            wr_p = dot(obj.fun_basis_Np(us),DM(obj.Pwr));
            tw_p = dot(obj.fun_basis_Np(us),DM(obj.Ptw));
            pos_sym = [x_p;y_p;z_p;wl_p;wr_p;tw_p];
            obj.fun_pos = Function('fun_pos',{us},{pos_sym},{'us'},{'pos_sym'});
        end
        function length_tot = length_eval(obj,uu)
            import casadi.*
            pointl = obj.fun_pos(uu);
            length_tot = sum(sqrt(diff(pointl(1,:)).^2 + diff(pointl(2,:)).^2 + diff(pointl(3,:)).^2));
            %disp(length_tot)
        end
        %% evaluate quadratic error
        function residual(obj)
            % method to evaluate the residual quadratic error
            import casadi.*
            if isempty(obj.Px) 
                error('First run obj.solve')
            end
            if isempty(obj.fun_pos)
                error('First run obj.eval')
            end
            point_middle = obj.fun_pos(obj.ub)';
            quadratic_error = (obj.x - point_middle(:,1)).^2 + (obj.y - point_middle(:,2)).^2 + (obj.z - point_middle(:,3)).^2;
            total_error = sum(quadratic_error);
            obj.q_error = quadratic_error;
            obj.tot_error = total_error;
            %disp(full(total_error));          
        end
        %% Basis for first order derivative
        function dbasis(obj)
            % build the basis function for NURBS derivative
            import casadi.*
            obj.Nd = SX(obj.n+2,obj.p+1);
            obj.Nd(:,1) = DM(1,1);
            if isempty(obj.N)
                error('first run obj.eval_symbolic()')
            end
            for j = 1:obj.p    
                for i = 1:(obj.n+1)  
                    num1 = if_else(DM(obj.u(i+j) - obj.u(i)) == 0,DM(1,1),j);
                    den1 = if_else(DM(obj.u(i+j) - obj.u(i)) == 0,1, obj.u(i+j) - obj.u(i));
                    num2 = if_else(DM(obj.u(i+j+1) - obj.u(i+1)) == 0,DM(1,1),j);
                    den2 = if_else(DM(obj.u(i+j+1) - obj.u(i+1)) == 0,1,obj.u(i+j+1) - obj.u(i+1));
                    obj.Nd(i,j+1) = obj.N(i,j)*(num1/den1)  - obj.N(i+1,j)*(num2/den2);
                end      
            end
            obj.Nd(obj.n+2,:) = DM(1,1);
            obj.Ndp = obj.Nd(1:obj.n+1,obj.p+1);
        end
        %% evaluate tangent vector (S-frame)
        function tangent(obj,us)
            % build the casadi function fun_ts(us) necessary to evaluate spline tangent vector
            import casadi.*
            if isempty(obj.Px)
                error('first run obj.solve')
            end
            if isempty(obj.fun_basis_Np)
                obj.basis(us);
            end
            obj.dbasis
            tx = dot(obj.Ndp,obj.Px);
            ty = dot(obj.Ndp,obj.Py);
            tz = dot(obj.Ndp,obj.Pz);
            obj.ts = [tx;ty;tz]./sqrt(tx^2+ty^2+tz^2);
            obj.fun_ts = Function('fun_ts',{us},{obj.ts},{'us'},{'ts'});
            norm_ts = sqrt(tx^2+ty^2+tz^2);
            obj.fun_normts = Function('fun_normts',{us},{norm_ts},{'us'},{'norm_ts'});
        end
        %% Basis for second order derivative
        function ddbasis(obj)
            % build basis function for NURBS second order derivative
            import casadi.*
            obj.Ndd = SX(obj.n+2,obj.p+1);
            obj.Ndd(:,1) = DM(1,1);
            if isempty(obj.Nd)
                error('first run obj.dbasis()')
            end
            for j = 1:obj.p    
                for i = 1:(obj.n+1)  
                    num1 = if_else(DM(obj.u(i+j) - obj.u(i)) == 0,DM(1,1),j);
                    den1 = if_else(DM(obj.u(i+j) - obj.u(i)) == 0,1, obj.u(i+j) - obj.u(i));
                    num2 = if_else(DM(obj.u(i+j+1) - obj.u(i+1)) == 0,DM(1,1),j);
                    den2 = if_else(DM(obj.u(i+j+1) - obj.u(i+1)) == 0,1,obj.u(i+j+1) - obj.u(i+1));
                    obj.Ndd(i,j+1) = obj.Nd(i,j)*(num1/den1) - obj.Nd(i+1,j)*(num2/den2);
                end      
            end
            obj.Ndd(obj.n+2,:) = DM(1,1);
            obj.Nddp = obj.Ndd(1:obj.n+1,obj.p+1);
            
        end
        %% evaluate normal vector (n of frenet-frame)
        function normal(obj,us)
            % build two casadi function, fun_kv(us) necessary to evaluate curvature without sign and fun_nv(us) necessary to evaluate theoretical normale vector
            import casadi.*
            if isempty(obj.Px)
                error('first run obj.solve')
            end
            if isempty(obj.fun_basis_Np)
                obj.basis(us);
            end
            obj.dbasis
            if isempty(obj.fun_ts)
                obj.tangent(us);
            end
            obj.ddbasis
            obj.nv = SX.sym('nv',3,1);
            nxv = dot(obj.Nddp,obj.Px);
            nyv = dot(obj.Nddp,obj.Py);
            nzv = dot(obj.Nddp,obj.Pz);
            obj.kv = [nxv;nyv;nzv]./(obj.fun_normts(us))^2 - (obj.fun_ts(us)./(obj.fun_normts(us))^2)*(dot(obj.fun_ts(us),[nxv;nyv;nzv]));
            obj.k = sqrt(obj.kv(1)^2 + obj.kv(2)^2 + obj.kv(3)^2);
            obj.nv = obj.kv./obj.k;
            obj.fun_kv = Function('fun_kv',{us},{obj.kv},{'us'},{'kv'});
            obj.fun_k  = Function('fun_k',{us},{obj.k},{'us'},{'k'});
            obj.fun_nv = Function('fun_nv',{us},{obj.nv},{'us'},{'nv'});
        end
        %% evaluate rotate normal vector (v of H-frame)
        function rnormal(obj,us)
            % build two casadi function fun_kp and fun_pt necessary to evaluate curvature with sign and rotate normal vector (vector normal to t in the x-y plane)
            import casadi.*
            if isempty(obj.Px)
                error('first run obj.solve');
            end
            if isempty(obj.fun_basis_Np)
                obj.basis(us);
            end
            obj.dbasis				 				   			 
            if isempty(obj.fun_ts)
                obj.tangent(us);
            end
            obj.ddbasis
            tvec = obj.fun_ts(us);
            vtx = -tvec(2);
            vty = tvec(1);
            vtz = 0*tvec(3);
            obj.vh = [vtx;vty;vtz]/sqrt(vtx^2 + vty^2);
            nxv = dot(obj.Nddp,obj.Px);
            nyv = dot(obj.Nddp,obj.Py);
            nzv = dot(obj.Nddp,obj.Pz);
            normal = [nxv;nyv;nzv]./(obj.fun_normts(us))^2 - (obj.fun_ts(us)./(obj.fun_normts(us))^2)*(dot(obj.fun_ts(us),[nxv;nyv;nzv]));
            obj.kp = dot(normal,obj.vh);
            obj.fun_kp  = Function('fun_kp',{us},{obj.kp},{'us'},{'kp'});
            obj.fun_vh = Function('fun_vh',{us},{obj.vh},{'us'},{'vh'});
        end
        %% evaluate binormal vector (b of frenet-frame)
        function binormal(obj,us)
            % method to build casadi function fun_bnv necessary to evaluate binormal vector from tangent and theoretical normal 
            import casadi.*
            if isempty (obj.fun_ts)
                obj.tangent(us)
            end
            if isempty(obj.fun_nv)
                obj.normal(us)
            end
            obj.bnv = cross(obj.fun_ts(us),obj.fun_nv(us));
            obj.fun_bnv = Function('fun_bnv',{us},{obj.bnv},{'us'},{'bnv'});
        end
        %% evaluate rbinormal vector (w of H-frame)
        function rbinormal(obj,us)
            % method to build casadi function fun_rbnv necessary to evaluate binormal vector from tangent and rotate normal 
            import casadi.*
            if isempty (obj.fun_ts)
                obj.tangent(us)
            end
            if isempty(obj.fun_vh)
                obj.rnormal(us)
            end
            obj.wh = cross(obj.fun_ts(us),obj.fun_vh(us));
            obj.fun_wh = Function('fun_wh',{us},{obj.wh},{'us'},{'wh'});
        end
		%% evaluate matrix Rgs-Rgh-Rhs & unit vector (ns, ms of S-frame) & g_gs & curvature along ns & geometric twist 
        function Sframe(obj,us)
            % build casadi function for Rgs-Rgh-Rhs-ns-ms-g_gs-kns-Vgs_s
            import casadi.*
            % run necessary method
            %obj.basis(us)
            obj.eval(us)
            obj.dbasis
            obj.ddbasis
            obj.tangent(us)
            obj.rnormal(us)
            obj.rbinormal(us)
            % evaluate position and H-frame vector
            pos = obj.fun_pos(us);
            eta = pos(6);
            t = obj.fun_ts(us);
            v = obj.fun_vh(us);
            w = obj.fun_wh(us);
            normts = obj.fun_normts(us);
            % evaluate banking
            %mu = atan(tan(eta)*cos(atan(-t(3)/sqrt(t(2)^2 + t(1)^2))));
            mu = eta;
            % Rotation matrix
            cmu = cos(mu);
            smu = sin(mu);
            Rhs = [1,0,0;0,cmu,-smu;0,smu,cmu];
            Rgh = [t(1),v(1),w(1);t(2),v(2),w(2);t(3),v(3),w(3)];
            obj.fun_Rgh = Function('fun_Rgh',{us},{Rgh},{'us'},{'Rgh'});
            Rgs = [t(1),v(1),w(1);t(2),v(2),w(2);t(3),v(3),w(3)]*[1,0,0;0,cmu,-smu;0,smu,cmu];
            obj.fun_Rgs = Function('fun_Rgs',{us},{Rgs},{'us'},{'Rgs'});
            % Vector of S-frame
            obj.fun_ns = Function('fun_ns',{us},{Rgs(:,2)},{'us'},{'ns'});
            obj.fun_ms = Function('fun_ms',{us},{Rgs(:,3)},{'us'},{'ms'});
            % Transformation matrix from G to S
            g_gs = [Rgs,[pos(1);pos(2);pos(3)];0,0,0,1];
            obj.fun_g_gs = Function('fun_g_gs',{us},{g_gs},{'us'},{'g_gs'});
            % Curvature in ns direction
            nsxv = dot(obj.Nddp,obj.Px);
            nsyv = dot(obj.Nddp,obj.Py);
            nszv = dot(obj.Nddp,obj.Pz);
            nts_squared = normts^2;
            nsnormal = [nsxv;nsyv;nszv]./nts_squared - (t./nts_squared)*(dot(t,[nsxv;nsyv;nszv]));            
            obj.kns = dot(nsnormal,Rgs(:,2));
            obj.fun_kns  = Function('fun_kns',{us},{obj.kns},{'us'},{'kns'});
            % Omega g-h^h
            b = dot(nsnormal,v);
            omega_x = b*t(3)/sqrt(t(1)^2 + t(2)^2);
            c = dot(nsnormal,w);
            omega_y = -c;
            omega_z = b;
            omega_gh_h = [0, - omega_z, omega_y; omega_z, 0, -omega_x; -omega_y, omega_x, 0];
            % omega h-s^h
            mu_dot = dot(obj.Ndp,obj.Ptw)/normts;
            obj.fun_mudot = Function('fun_mudot',{us},{mu_dot},{'us'},{'mu_dot'});
            omega_hs_h = [0, 0, 0;0, 0, -mu_dot;0, mu_dot, 0];
            % omega g-s^s
            omega_gs_s = Rhs'*(omega_gh_h + omega_hs_h)*Rhs;
            omega_gs_sv = [omega_gs_s(3,2);omega_gs_s(1,3);omega_gs_s(2,1)];
            % tangent vector in S-frame
            ts_s = Rgs'*t;
            % Twist vector in S-frame
            Vgs_s = [ts_s; omega_gs_sv];
            obj.fun_Vgs_s = Function('fun_Vgs_s',{us},{Vgs_s},{'us'},{'Vgs_s'});
        end		
        %% plot track_point vs nurbs
        function plot_track_vs_nurbs(obj)
            % method to plot track and nurbs evaluated in the ub points
            import casadi.*
            obj.solve
            obj.eval(SX.sym('us'))
            obj.tangent(SX.sym('us'))
            obj.normal(SX.sym('us'))
            obj.binormal(SX.sym('us'))
            pos = full(obj.fun_pos(obj.ub));
            plot3(pos(1,:),pos(2,:),pos(3,:))
            hold on
            plot3(pos(1,:),pos(2,:),pos(3,:),'o')
            hold on
            plot3(obj.x, obj.y, obj.z, 'x')
        end
        %% plot with tangent and norm
        function plot_tn(obj,us)
            % method to plot 3D middle line with tangent, normal and binormal versor (from curve theory)
            import casadi.*
            obj.solve
            obj.eval(SX.sym('us'))
            obj.tangent(SX.sym('us'))
            obj.normal(SX.sym('us'))
            obj.binormal(SX.sym('us'))
            pos = full(obj.fun_pos(us));
            tg = full(obj.fun_ts(us));
            ng = full(obj.fun_nv(us));
            bng = full(obj.fun_bnv(us));
            quiver3(pos(1,:),pos(2,:),pos(3,:),tg(1,:),tg(2,:),tg(3,:))
            axis equal
            hold on
            quiver3(pos(1,:),pos(2,:),pos(3,:),ng(1,:),ng(2,:),ng(3,:),0.5)
            axis equal
            hold on
            quiver3(pos(1,:),pos(2,:),pos(3,:),bng(1,:),bng(2,:),bng(3,:),0.5)
            axis equal
            hold on
            plot3(pos(1,:),pos(2,:),pos(3,:))
        end
        %% plot with tangent and rnorm
        function plot_trn(obj,us)
            % method to plot 3D middle line with tangent, rotate normal and rotate binormal versor
            import casadi.*
            obj.solve
            obj.eval(SX.sym('us'))
            obj.Sframe(SX.sym('us'));
            pos = full(obj.fun_pos(us));
            tg = full(obj.fun_ts(us));
            vg = full(obj.fun_vh(us));
            wg = full(obj.fun_wh(us));
            ng = full(obj.fun_ns(us));
            mg = full(obj.fun_ms(us));
            quiver3(pos(1,:),pos(2,:),pos(3,:),tg(1,:),tg(2,:),tg(3,:))
            axis equal
            hold on
            quiver3(pos(1,:),pos(2,:),pos(3,:),vg(1,:),vg(2,:),vg(3,:),0.5)
            axis equal
            hold on
            quiver3(pos(1,:),pos(2,:),pos(3,:),wg(1,:),wg(2,:),wg(3,:),0.5)
            axis equal
            hold on
            quiver3(pos(1,:),pos(2,:),pos(3,:),ng(1,:),ng(2,:),ng(3,:),0.5)
            axis equal
            hold on
            quiver3(pos(1,:),pos(2,:),pos(3,:),mg(1,:),mg(2,:),mg(3,:),0.5)
            axis equal
            hold on
            plot3(pos(1,:),pos(2,:),pos(3,:))
        end
        %% complete evaluation tangent and norm
        function full_eval_norm(obj,us)
            % method to evaluate all method to build NURBS (t, n, bnv)
            import casadi.*
            obj.solve
            obj.eval(us)
            obj.tangent(us)
            obj.normal(us)
            obj.binormal(us)
        end
        %% complete evaluation tangent and rnorm
        function full_eval_rnorm(obj,us)
            % method to evaluate all method to build NURBS (t, vh, wh)
            import casadi.*
            obj.solve
            obj.eval(us)
            %obj.tangent(us)
            %obj.rnormal(us)
            %obj.rbinormal(us)
			obj.Sframe(us)	 
        end
        %% Plot 3D track with boundaries
        function ax = full_plot(obj, us, options)
            arguments 
                obj;
                us;
                options.parent = [];
                options.patch = true;
                options.dotted = false;
                options.closed = false;
            end
            parent = options.parent;
            % method to plot 3D track with boundaries
            import casadi.*
            if isempty(obj.fun_Rgs)
                obj.full_eval_rnorm(SX.sym('us'));
            end
            Point = zeros(length(us), 9);
            middle = obj.fun_pos(us);
            xm = full(middle(1, :));
            ym = full(middle(2, :));
            zm = full(middle(3, :));
            wm = full(middle(4, :));
            mum = middle(6);
            tm = obj.fun_ts(us);
            Rgs = reshape(full(obj.fun_Rgs(us)), 3, 3, length(us));
            ext = zeros(3, length(us));
            int = zeros(3, length(us));
            for i=1:length(us)
                ext(:, i) = full(Rgs(:, 2, i).*wm(i));
                int(:, i) = full(Rgs(:, 2, i).*-wm(i));       
            end
            Point(:, 4) = xm ;
            Point(:, 5) = ym ;
            Point(:, 6) = zm;
            Point(:,1:3) = ([xm; ym; zm+3e-3] + ext)';
            Point(:,7:9) = ([xm; ym; zm+3e-3] + int)';
            
            obj.Points = Point;
            if isempty(parent)
                figure; axis('equal'); hold on; grid on
                parent = gca;
            end
            id_end = length(xm);
            if ~options.closed
                id_end = id_end-1;
            end
            if options.dotted
                plot3(Point(1:id_end,1),Point(1:id_end,2),Point(1:id_end,3)+3e-3,'r','parent', parent,'Linewidth',1., 'Linestyle','none','marker','*','HandleVisibility','off');
                plot3(Point(1:id_end,4),Point(1:id_end,5),Point(1:id_end,6)+3e-3,'g','parent', parent,'Linewidth',0.5, 'Linestyle','none','marker','*','HandleVisibility','off');
                plot3(Point(1:id_end,7),Point(1:id_end,8),Point(1:id_end,9)+3e-3,'r','parent', parent,'Linewidth',1., 'Linestyle','none','marker','*','HandleVisibility','off')
            else
                plot3(Point(1:id_end,1),Point(1:id_end,2),Point(1:id_end,3)+3e-3,'color', 'black','parent', parent,'Linewidth',1,'HandleVisibility','off');
                plot3(Point(1:id_end,4),Point(1:id_end,5),Point(1:id_end,6)+1,'color', 'black','parent', parent,'Linewidth',1, 'linestyle', '--','DisplayName','Middle Line');
                plot3(Point(1:id_end,7),Point(1:id_end,8),Point(1:id_end,9)+3e-3,'color', 'black','parent', parent,'Linewidth',1,'HandleVisibility','off');
            end
            
            zmin = max(abs(zm));
            
            if options.patch
                surface('XData', [Point(1:id_end,1),(Point(1:id_end,7))], 'YData', [Point(1:id_end,2),(Point(1:id_end,8))], 'ZData',[Point(1:id_end,3),(Point(1:id_end,9))], 'facecolor', [0.6, 0.6, 0.6], 'parent', parent, 'edgecolor', 'none','HandleVisibility','off')
                surface('XData', [Point(1:id_end,1),(Point(1:id_end,7))], 'YData', [Point(1:id_end,2),(Point(1:id_end,8))], 'ZData',[Point(1:id_end,3),(Point(1:id_end,9))].*0 - zmin, 'facecolor', [146 36 40]./255, 'parent', parent,'HandleVisibility','off')
                surface('XData', [Point(1:id_end,1),(Point(1:id_end,1))], 'YData', [Point(1:id_end,2),(Point(1:id_end,2))], 'ZData',[Point(1:id_end,3),Point(1:id_end,3).*0-zmin]-1e-4, 'facecolor', [146 36 40]./255, 'parent', parent, 'edgecolor', 'k','HandleVisibility','off')
                surface('XData', [Point(1:id_end,7),(Point(1:id_end,7))], 'YData', [Point(1:id_end,8),(Point(1:id_end,8))], 'ZData',[Point(1:id_end,9),Point(1:id_end,9).*0-zmin]-1e-4, 'facecolor', [146 36 40]./255, 'parent', parent, 'edgecolor', 'k','HandleVisibility','off')
            end
            if nargout>0
               ax = parent; 
            end
        end
        function ax = full_plot_admm(obj,us, options)
            arguments 
                obj;
                us;
                options.parent = [];
                options.patch = true;
                options.dotted = false;
                options.closed = false;
            end
            parent = options.parent;
            % method to plot 3D track with boundaries
            import casadi.*
            if isempty(obj.fun_Rgs)
                obj.full_eval_rnorm(SX.sym('us'));
            end
            Point = zeros(length(us), 9);
            middle = obj.fun_pos(us);
            xm = full(middle(1, :));
            ym = full(middle(2, :));
            zm = full(middle(3, :));
            wm = 2*full(middle(4, :));
            mum = middle(6);
            tm = obj.fun_ts(us);
            Rgs = reshape(full(obj.fun_Rgs(us)), 3, 3, length(us));
            ext = zeros(3, length(us));
            int = zeros(3, length(us));
            for i=1:length(us)
                ext(:, i) = full(Rgs(:, 2, i).*wm(i));
                int(:, i) = full(Rgs(:, 2, i).*-wm(i));       
            end
            Point(:, 4) = xm ;
            Point(:, 5) = ym ;
            Point(:, 6) = zm;
            Point(:,1:3) = ([xm; ym; zm+3e-3] + ext)';
            Point(:,7:9) = ([xm; ym; zm+3e-3] + int)';
            
            obj.Points = Point;
            if isempty(parent)
                figure; axis('equal'); hold on; grid on
                parent = gca;
            end
            id_end = length(xm);
            if ~options.closed
                id_end = id_end-1;
            end
            
            if options.patch
                X = [Point(1:id_end,1),(Point(1:id_end,7))];
                Y = [Point(1:id_end,2),(Point(1:id_end,8))];
                 fill( X(:), Y(:), 0*Y(:) , 'facecolor', [0.6, 0.6, 0.6], 'parent', parent, 'edgecolor', 'none','HandleVisibility','off')                
%                 surface('XData', [Point(1:id_end,1),(Point(1:id_end,7))], 'YData', [Point(1:id_end,2),(Point(1:id_end,8))], 'ZData',[Point(1:id_end,3),(Point(1:id_end,9))], 'facecolor', [0.6, 0.6, 0.6], 'parent', parent, 'edgecolor', 'none','HandleVisibility','off')
%                 surface('XData', [Point(1:id_end,1),(Point(1:id_end,7))], 'YData', [Point(1:id_end,2),(Point(1:id_end,8))], 'ZData',[Point(1:id_end,3),(Point(1:id_end,9))].*0 - zmin, 'facecolor', [146 36 40]./255, 'parent', parent,'HandleVisibility','off')
%                 surface('XData', [Point(1:id_end,1),(Point(1:id_end,1))], 'YData', [Point(1:id_end,2),(Point(1:id_end,2))], 'ZData',[Point(1:id_end,3),Point(1:id_end,3).*0-zmin]-1e-4, 'facecolor', [146 36 40]./255, 'parent', parent, 'edgecolor', 'k','HandleVisibility','off')
%                 surface('XData', [Point(1:id_end,7),(Point(1:id_end,7))], 'YData', [Point(1:id_end,8),(Point(1:id_end,8))], 'ZData',[Point(1:id_end,9),Point(1:id_end,9).*0-zmin]-1e-4, 'facecolor', [146 36 40]./255, 'parent', parent, 'edgecolor', 'k','HandleVisibility','off')
            end
            if options.dotted
                plot(Point(1:id_end,1),Point(1:id_end,2),'r','parent', parent,'Linewidth',1., 'Linestyle','none','marker','*','HandleVisibility','off');
                plot(Point(1:id_end,4),Point(1:id_end,5),'g','parent', parent,'Linewidth',0.5, 'Linestyle','none','marker','*','HandleVisibility','off');
                plot(Point(1:id_end,7),Point(1:id_end,8),'r','parent', parent,'Linewidth',1., 'Linestyle','none','marker','*','HandleVisibility','off')
            else
                plot(Point(1:id_end,1),Point(1:id_end,2),'color', 'black','parent', parent,'Linewidth',1,'HandleVisibility','off');
                plot(Point(1:id_end,4),Point(1:id_end,5),'color', 'black','parent', parent,'Linewidth',1, 'linestyle', '--','DisplayName','Middle Line');
                plot(Point(1:id_end,7),Point(1:id_end,8),'color', 'black','parent', parent,'Linewidth',1,'HandleVisibility','off');
            end
            
            zmin = max(abs(zm));
            

            if nargout>0
               ax = parent; 
            end
        end
            
        %% evaluate symbolic functions with smaller computational graph size
        % this is achieved knowing the numerical value of the nurbs
        % variable that the optimization should fall in.
        % The numerical value then is used to evaluate the only necessary
        % knot span neighbourhood
        
        function tangent_subGraph(obj)
            % build the casadi function fun_ts(us) necessary to evaluate spline tangent vector
            import casadi.*
            if isempty(obj.Px)
                error('first run obj.solve')
            end
            tx = dot(obj.Ndp,obj.Px);
            ty = dot(obj.Ndp,obj.Py);
            tz = dot(obj.Ndp,obj.Pz);
            norm_ts = sqrt(tx^2+ty^2+tz^2);
            obj.ts = [tx;ty;tz]./norm_ts;
            obj.nts = norm_ts;
%             obj.fun_ts = Function('fun_ts',{us},{obj.ts},{'us'},{'ts'});
%             obj.fun_normts = Function('fun_normts',{us},{norm_ts},{'us'},{'norm_ts'});
        end
        
        function  eval_subGraph(obj,us)
            % build the casadi function fun_pos(us) necessary to evaluate the spline
            % insert us = SX.sym('us')
            %if Barcellona is our instance we run
            %Barcellona.eval_symbolic(SX.sym('us')) and then we have
            %Barcellona.fun_pos(us) with us numerical or symbolic to obtain
            %symbolic or numerical position. us could be a vector.
            import casadi.*
            if isempty(obj.Px)
                error('First run obj.solve')
            end
            x_p = dot(obj.Np,obj.Px);
            y_p = dot(obj.Np,obj.Py);
            z_p = dot(obj.Np,obj.Pz);
            wl_p = dot(obj.Np,obj.Pwl);
            wr_p = dot(obj.Np,obj.Pwr);
            tw_p = dot(obj.Np,obj.Ptw);
            pos_sym = [x_p;y_p;z_p;wl_p;wr_p;tw_p];
            obj.fun_pos = Function('fun_pos',{us},{pos_sym},{'us'},{'pos_sym'});
        end
        
        function rnormal_subGraph(obj)
            % build two casadi function fun_kp and fun_pt necessary to evaluate curvature with sign and rotate normal vector (vector normal to t in the x-y plane)
            import casadi.*
            if isempty(obj.Px)
                error('first run obj.solve');
            end			 				   			 
            tvec = obj.ts;
            norm_ts = obj.nts;
            vtx = -tvec(2);
            vty = tvec(1);
            vtz = 0*tvec(3);
            obj.vh = [vtx;vty;vtz]/sqrt(vtx^2 + vty^2);
            nxv = dot(obj.Nddp,obj.Px);
            nyv = dot(obj.Nddp,obj.Py);
            nzv = dot(obj.Nddp,obj.Pz);
            normal = [nxv;nyv;nzv]./(norm_ts)^2 - (tvec./(norm_ts)^2)*(dot(tvec,[nxv;nyv;nzv]));
            obj.kp = dot(normal,obj.vh);
%             obj.fun_kp  = Function('fun_kp',{us},{obj.kp},{'us'},{'kp'});
%             obj.fun_vh = Function('fun_vh',{us},{obj.vh},{'us'},{'vh'});
        end
        
        function Sframe_subGraph(obj,us)
            % build casadi function for Rgs-Rgh-Rhs-ns-ms-g_gs-kns-Vgs_s
            import casadi.*
            % evaluate position and H-frame vector
            pos = obj.fun_pos(us);
            eta = pos(6);
%             t = obj.fun_ts(us);
%             v = obj.fun_vh(us);
%             w = obj.fun_wh(us);
            t = obj.ts;
            v = obj.vh;
            w = obj.wh;
            norm_ts = obj.nts;
            % evaluate banking
            %mu = atan(tan(eta)*cos(atan(-t(3)/sqrt(t(2)^2 + t(1)^2))));
            mu = eta;
            % Rotation matrix
            cmu = cos(mu);
            smu = sin(mu);
            Rhs = [1,0,0;0,cmu,-smu;0,smu,cmu];

            Rgs = [t(1),v(1),w(1);t(2),v(2),w(2);t(3),v(3),w(3)]*[1,0,0;0,cmu,-smu;0,smu,cmu];
            % Vector of S-frame

            % Transformation matrix from G to S
            g_gs = [Rgs,[pos(1);pos(2);pos(3)];0,0,0,1];
            obj.fun_g_gs = Function('fun_g_gs',{us},{g_gs},{'us'},{'g_gs'});
            % Curvature in ns direction
            nsxv = dot(obj.Nddp,obj.Px);
            nsyv = dot(obj.Nddp,obj.Py);
            nszv = dot(obj.Nddp,obj.Pz);
            nsnormal = [nsxv;nsyv;nszv]./(obj.fun_normts(us))^2 - (obj.fun_ts(us)./(obj.fun_normts(us))^2)*(dot(obj.fun_ts(us),[nsxv;nsyv;nszv]));            
            obj.kns = dot(nsnormal,Rgs(:,2));
%             obj.fun_kns  = Function('fun_kns',{us},{obj.kns},{'us'},{'kns'});
            % Omega g-h^h
            b = dot(nsnormal,v);
            omega_x = b*sqrt(1- t(1)^2 - t(2)^2)/sqrt(t(1)^2 + t(2)^2);
            c = dot(nsnormal,w);
            omega_y = -c;
            omega_z = b;
            omega_gh_h = [0, - omega_z, omega_y; omega_z, 0, -omega_x; -omega_y, omega_x, 0];
            % omega h-s^h
            mu_dot = dot(obj.Ndp,obj.Ptw)/norm_ts;
%             obj.fun_mudot = Function('fun_mudot',{us},{mu_dot},{'us'},{'mu_dot'});
            omega_hs_h = [0, 0, 0;0, 0, -mu_dot;0, mu_dot, 0];
            % omega g-s^s
            omega_gs_s = Rhs'*(omega_gh_h + omega_hs_h)*Rhs;
            omega_gs_sv = [omega_gs_s(3,2);omega_gs_s(1,3);omega_gs_s(2,1)];
            % tangent vector in S-frame
            ts_s = Rgs'*t;
            % Twist vector in S-frame
            Vgs_s = [ts_s; omega_gs_sv];
            obj.fun_Vgs_s = Function('fun_Vgs_s',{us},{Vgs_s},{'us'},{'Vgs_s'});
        end		
        
        function rbinormal_subGraph(obj)
            % method to build casadi function fun_rbnv necessary to evaluate binormal vector from tangent and rotate normal 

            obj.wh = cross(obj.ts,obj.vh);
%             obj.fun_wh = Function('fun_wh',{us},{obj.wh},{'us'},{'wh'});
        end
        
        function update_subGraphBases(obj, us, us_num, deltaU)
            %
            % deltaU = integration step size of the nurbs variable. it is
            % needed to robustly detect if neighbours knot spans should be
            % added for numerical stability during optimization
            %
            import casadi.*
            % this function choses only the neighbour knot spans based on
            % the numerical value that the function should be close to...
            basisDerivatives = subGraph_DersBasisFuns(us_num, obj.p, obj.u, 2, deltaU);
            N_der = basisDerivatives(us)';
            obj.Np = N_der(:, 1);
            obj.Ndp = N_der(:, 2);
            obj.Nddp = N_der(:, 3);
            obj.tangent_subGraph();
            obj.rnormal_subGraph();
            obj.rbinormal_subGraph();
            eval_subGraph(obj,us);
            Sframe_subGraph(obj,us);
        end
        
    end
           
end
        

            

                    




  
