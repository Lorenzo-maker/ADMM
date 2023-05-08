%% admm settings
% The settings written here are used by the ADMM instance and the sub
% problems instances

%%%%%%%%%%%%%%%%%%%%% Sub-problem parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lap = 1;            % number of lap
Nproblems = 4*lap;  % number of subproblems per lap
Nsteps = 1500*lap;  % TOTAL STEPS 
e = 30;             % if e > 0 it must be greater than o/2 % number of mesh interval for overlapping area
o = 1;              % number of mesh interval for consensus (if o = 0 consensus only on states at interfaces)
d = 2;              % number of collocation points 

%%%%%% check discretization type and define final alpha %%%%%%       
alpha_numeric = true;

global alfa_end
alfa_end = 1;
dalfa = alfa_end/(Nsteps/lap);

%%%%%%%% Scaling for variables and variables size %%%%%%%%%%%%%%%

car_parameters_ocp; 

%%%%%%%%%%%%%%%%%%%%% rho parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho_scale_head = 1;%0.01; % scaling factor for rho initialization   
rho_scale_tail = 1;%0.02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% Initial guess %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
init_guess = true;
if init_guess 
    load guess_0_1_tot_tris.mat guess_0_1_tot_tris % guess_0_1colloc.mat for whole track
    guess = guess_0_1_tot_tris;    
end

%%%%%%%%%%%%%%%%%%%%% Update Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
update_method = 'ADMM'; %ADMM, Stationary, Nesterov, Automatic
alfa_stationary = 0.3;
eta = 0.85;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileparts(which('guess_0_1_tot_tris.mat'))
%%%%%%%%%%%%%%%%%%%%% Tolerance for X, U, Z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon_state = 0.01*ones(nx,1);      % define error tolerance for states
epsilon_control = 0.001*ones(nu,1);   % define error tolerance for controls
epsilon_algebraic = 0.001*ones(nz,1); % define error tolerance for algebraic parameters
error_factor = 0.5;                   % coefficient to scale epsilon to obtain tolerance for error |X-Z|
delta_z_factor = 0.25;                % coefficient to scale epsilon to obtain tolerance for |Zk+1 - Zk|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% create tollerance array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if o > 0 
    vec = [epsilon_control;epsilon_algebraic;epsilon_state];
    epsilon_array = [epsilon_state;repmat(vec,o,1)];
else
    epsilon_array = epsilon_state;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Nproblems == 1 % single problem
    % IPOPT options
    IPOPT_opt = struct;
    IPOPT_opt.ipopt.linear_solver = 'ma57';
    IPOPT_opt.ipopt.ma57_pre_alloc = 10;
    IPOPT_opt.ipopt.max_iter = 800;
    IPOPT_opt.ipopt.print_level = 5;
    IPOPT_opt.print_time = 1;
    IPOPT_opt.ipopt.mu_init = 1e-3;
%   IPOPT_opt.ipopt.warm_start_init_point = 'yes';
%   IPOPT_opt.ipopt.alpha_for_y = 'min-dual-infeas';
%   IPOPT_opt.ipopt.nlp_scaling_method = 'none';
    IPOPT_opt.ipopt.tol = 1e-6;
    IPOPT_opt.expand = false;
%   IPOPT_opt.ipopt.hessian_approximation = 'limited-memory';
%   IPOPT_opt.ipopt.alpha_min_frac = 0.6;
%   IPOPT_opt.ipopt.nlp_scaling_max_gradient = 10^4;
%   IPOPT_opt.ipopt.mu_strategy = 'adaptive';
%   IPOPT_opt.ipopt.mu_oracle = 'probing';
%   IPOPT_opt.ipopt.adaptive_mu_globalization = 'never-monotone-mode';
else % ADMM
    % IPOPT options
    IPOPT_opt = struct;
    IPOPT_opt.ipopt.linear_solver = 'ma57';
    IPOPT_opt.ipopt.ma57_pre_alloc = 10;
    IPOPT_opt.ipopt.max_iter = 800;
    IPOPT_opt.ipopt.print_level = 5;
    IPOPT_opt.print_time = 1;
    IPOPT_opt.ipopt.mu_init = 1e-3;
    IPOPT_opt.ipopt.warm_start_init_point = 'yes';
%   IPOPT_opt.ipopt.nlp_scaling_method =  'none';
%   IPOPT_opt.ipopt.alpha_for_y = 'min-dual-infeas';
%   IPOPT_opt.ipopt.nlp_scaling_method = 'none';
    IPOPT_opt.ipopt.tol = 1e-6;
    IPOPT_opt.expand = false;
%   IPOPT_opt.ipopt.alpha_min_frac = 0.6;
%   IPOPT_opt.ipopt.nlp_scaling_max_gradient = 10^4;
%   IPOPT_opt.ipopt.mu_strategy = 'adaptive';
    IPOPT_opt.ipopt.ma57_small_pivot_flag = 1;
    IPOPT_opt.ipopt.ma57_pivtol = 1e-5;
%   IPOPT_opt.ipopt.mu_oracle = 'probing';
%   IPOPT_opt.ipopt.adaptive_mu_globalization = 'never-monotone-mode';
end