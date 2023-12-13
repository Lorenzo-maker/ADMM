function w0 = warm_start(self,track,d,N, index_start, index_end)

    % Initialize the MLTP solution. First compute the initial state and then propagate conditions to whole solution.
    %
    % NOTE Numeric values in this function are based on a FSAE car with the default parameters. Different car models may require some adaptation.

    % Set target velocity
    speed_target = 21;

    % Compute torque balancing the aerodynamic drag
    T_target = 0.5*self.data.aero_S*self.data.aero_rho*self.data.aero_Cx*speed_target^2*self.data.tire_radius{2};

    % Define the time base for ws simulaiton
    t_start = 0;
    t_push = t_start+1;
    t_release = t_push+speed_target/7; % NOTE Formula empirica raffinatissima.
    t_stop = t_release+1;
    h = 0.001;
    t = t_start:h:t_stop;
    M = length(t)-1;

    % Set the initial condition for warm-start simulation
    V_gb_init = [0;0;0;0;0;0];
    dt_trim_init = [0;0;0;0];
    omega_init = [0;0;0;0];
    q_gb_init = [0;0;0.25;0;0.0048;0];
    trim_init = [-0.003;-0.003;0.003;0.003];
    x_unscaled = [V_gb_init;dt_trim_init;omega_init;q_gb_init;trim_init];
    x_init = full(self.data.X_scale.\x_unscaled);

    % Design the (normalized) input trajectory stabilizing vehicle on target conditions
    T_unit = zeros(1,M+1);
    T_unit(t>t_push) = (t(t>t_push)-t_push)/0.25;
    T_unit(t>t_push+0.25) = 1;
    T_unit(t>t_release) = 1-(t(t>t_release)-t_release)/0.25;
    T_unit(t>t_release+0.25) = T_target/self.data.U_scale(1);
    u = [T_unit(2:end);zeros(1,M);zeros(1,M)];

    % Run simulation to find the steady-state conditions
    [x,v] = self.simulate(t_start,t_stop,h,x_init,u);

    % Plot simulation results to check that ss conditions are reached
    % plot(axes(figure('Name','States')),t,self.data.x_scale.*x)
    % plot(axes(figure('Name','Sparsifying variables')),t(2:end),v)

    % Filter out null ss components and mirror symmetric ones
    x_ss = x(:,end);
    x_ss([2,4,5,6,7,8,9,10,15,16,18,20]) = 0;
    x_ss([11,12]) = (x_ss(11)+x_ss(12))/2;
    x_ss([13,14]) = (x_ss(13)+x_ss(14))/2;
    x_ss([21,22]) = (x_ss(21)+x_ss(22))/2;
    x_ss([23,24]) = (x_ss(23)+x_ss(24))/2;
    v_ss = v(:,end);
    v_ss([1,2]) = (v_ss(1)+v_ss(2))/2;
    v_ss([3,4]) = (v_ss(3)+v_ss(4))/2;
    v_ss(5) = (v_ss(5)+v_ss(6))/2;
    v_ss(6) = -v_ss(5);
    v_ss(7) = (v_ss(7)+v_ss(8))/2;
    v_ss(8) = -v_ss(7);
    v_ss([9,10]) = (v_ss(9)+v_ss(10))/2;
    v_ss([11,12]) = (v_ss(11)+v_ss(12))/2;

    % Compute initial state given the track geometry and the ss conditions
    speed_ss = sqrt((x_ss(1)*self.data.X_scale(1))^2+(x_ss(3)*self.data.X_scale(3))^2);
    height_ss = self.data.X_scale(17)*x_ss(17);
    pitch_ss = self.data.X_scale(19)*x_ss(19);
    x0 = x_ss;
    R = reshape(self.data.track.d_gsR_gs(4:12,index_start),3,3)*rotY(pitch_ss);
    x0(15:17) = self.data.track.d_gsR_gs(1:3,index_start)+height_ss*R(:,3);
    x0(15:17) = self.data.X_scale(15:17).\x0(15:17);
    x0(18:20) = [atan2(R(2,1),R(1,1));atan2(-R(3,1),sqrt(1-R(3,1)^2));atan2(R(3,2),R(3,3))];
    x0(18) = self.data.track.psi_grid(index_start);
    tmp = x0(18);
    x0(18:20) = self.data.X_scale(18:20).\x0(18:20);

    % Extend ss conditions to whole solution
    w0 = nan(self.nu+self.nz+self.nx*(d+1),N);
    u_k = self.data.U_scale.\[T_target;0;0];
    points = self.data.track.d_gsR_gs(1:3,index_start:index_end);
    track_len = sum(sqrt(sum(diff(points')'.^2)));
    z_k = self.data.Z_scale(1:3).\[0;height_ss;track_len/N/speed_ss];
    z_k = [z_k;v_ss];
    x_k = x_ss;
    for k = 1:N
        R = reshape(self.data.track.d_gsR_gs(4:12,index_start + k),3,3)*rotY(pitch_ss);
        x_k(15:17) = self.data.track.d_gsR_gs(1:3,index_start + k)+height_ss*R(:,3);
        x_k(15:17) = self.data.X_scale(15:17).\x_k(15:17);
        x_k(18:20) = [atan2(R(2,1),R(1,1));atan2(-R(3,1),sqrt(1-R(3,1)^2));atan2(R(3,2),R(3,3))];
        x_k(18) = self.data.track.psi_grid(index_start+k);
%         e = x_k(18)-tmp;
%         if e>pi
%             x_k(18) = x_k(18)-2*pi*round(e/(2*pi));
%         elseif e<-pi
%             x_k(18) = x_k(18)+2*pi*round(-e/(2*pi));
%         end
        x_k(6) = (x_k(18)-tmp)/(self.data.Z_scale(3)*z_k(3));
        tmp = x_k(18);
        x_k(18:20) = self.data.X_scale(18:20).\x_k(18:20);
        x_k(6) = self.data.X_scale(6)\x_k(6);
        w_k = [u_k;z_k;repmat(x_k,d+1,1)];
        w0(:,k) = w_k;
    end
    w0 = reshape(w0,[],1);
    w0 = [x0;w0];

end