%% admm settings
% The settings written here are used by the ADMM instance and the sub
% problems instances

%%%%%%%%%%%%%%%%%%%%% Build Track %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%% Sub-problem parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lap = 1;
Nproblems = 4*lap;
Nsteps = 1500*lap;% TOTAL STEPS 
overlap = 81;%81;%81; % ***!!! dispari
overlap_control = true;
direct_c = true;
d = 2; % this collocation order is actually inside direct_collocation.m
       % changing it here does not take any effect (used only for state
       % size calculations)
global alfa_end
alfa_end = 1;%0.99;%0.99;
dalfa = alfa_end/(Nsteps/lap);

car_parameters_ocp;

%%%%%%%%%%%%%%%%%%%%%% variable sizes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = 11;
nu = 7;
nz = 8;

%%%%%%%%%%%%%%%%%%%%% rho parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mod_iter = 2; % rho increase/decrease iteration
rho_increase = 1.1; % rho increasing factor
rho_decrease = 0.95; % rho decreasing factor
rho_scale_head = 1;%0.01; % scaling factor for rho initialization   
rho_scale_tail = 1;%0.02;
rho_max_factor = 5; %maximum increase factor for rho respect to initial value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load guess_0_1_tot_bis.mat guess_0_1_tot_bis % guess_0_1colloc.mat for whole track
% guess = guess_0_1_tot_bis;
load guess_0_1_tot_tris.mat guess_0_1_tot_tris % guess_0_1colloc.mat for whole track
guess = guess_0_1_tot_tris;
% load guess_0_1_curvPID_v2.mat guess_0_1_curvPID % guess_0_1colloc.mat for whole track
% guess = guess_0_1_curvPID;
% guess.q = guess.q';
% guess.q = guess.q(:);
%%%%%%%%%%%%%%%%%%%%% Update Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
update_method = 'ADMM'; %ADMM, Stationary, Nesterov, Automatic
alfa_stationary = 0.3;
eta = 0.85;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% Tolerance for X, U, Z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon_state = 0.01*ones(nx,1);
epsilon_control = 0.001*ones(nu,1);
epsilon_algebraic = 0.001*ones(nz,1);
error_factor = 0.5;  % what is this ?
delta_z_factor = 0.25; % what is this ?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% create tollerance array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% only one consensus
elemOverlap = (1*(nx + nu + nz)+nx);

epsilon_array = [];
if elemOverlap > nx 
    while length(epsilon_array) < elemOverlap
        epsilon_array = [epsilon_array;epsilon_state;epsilon_control;epsilon_algebraic;epsilon_state];
    end
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
%         IPOPT_opt.ipopt.warm_start_init_point = 'yes';
    %     IPOPT_opt.ipopt.alpha_for_y = 'min-dual-infeas';
    %     IPOPT_opt.ipopt.nlp_scaling_method = 'none';
    IPOPT_opt.ipopt.tol = 1e-6;
    IPOPT_opt.expand = false;
%     IPOPT_opt.ipopt.hessian_approximation = 'limited-memory';
    %     IPOPT_opt.ipopt.alpha_min_frac = 0.6;
    %     IPOPT_opt.ipopt.nlp_scaling_max_gradient = 10^4;
%     IPOPT_opt.ipopt.mu_strategy = 'adaptive';
%     IPOPT_opt.ipopt.mu_oracle = 'probing';
%     IPOPT_opt.ipopt.adaptive_mu_globalization = 'never-monotone-mode';
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
    %IPOPT_opt.ipopt.nlp_scaling_method =  'none';

    %     IPOPT_opt.ipopt.alpha_for_y = 'min-dual-infeas';
    %     IPOPT_opt.ipopt.nlp_scaling_method = 'none';
    IPOPT_opt.ipopt.tol = 1e-6;
    IPOPT_opt.expand = false;
    %     IPOPT_opt.ipopt.alpha_min_frac = 0.6;
    %     IPOPT_opt.ipopt.nlp_scaling_max_gradient = 10^4;
%     IPOPT_opt.ipopt.mu_strategy = 'adaptive';
    IPOPT_opt.ipopt.ma57_small_pivot_flag = 1;
    IPOPT_opt.ipopt.ma57_pivtol = 1e-5;
    % IPOPT_opt.ipopt.mu_oracle = 'probing';
    % IPOPT_opt.ipopt.adaptive_mu_globalization = 'never-monotone-mode';
end