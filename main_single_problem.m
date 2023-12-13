%% ADMM main instance script
restoredefaultpath
clear all
clc
close all
working_dir = pwd;

%savepath = 'Results\1500_1_0'; % for results path

%% path

addpath(genpath('../casadi_3_6_4'));
addpath(genpath('Classes'));
addpath(genpath('Track'));
addpath(genpath('ADMM_functions'));
addpath(genpath('Data'));
addpath(genpath('scripts'));
addpath(genpath('utils'));
addpath(genpath('Extra'));

%% Create track
import casadi.*

%% Create track
% 1 -> load saved track otherwise build track
track_saved = 1;
%1 -> Catalunya, 2 -> Nurbrugring, 3 -> Calabogie
track_index = 2; 
% Build track
warning off
switch track_saved 
    case 1
        load('track_stored');
        track = track_stored{track_index};
    otherwise
        cntrlpts = 700; % Number of control points
        track = track_fun(track_index, '3D', 1, 'end', contrlpts);
end
pista = track;
pista.residual;
disp(['fitting error is ', num2str(full(pista.tot_error))])

save('Data\pista.mat', 'pista');
%% settings
ADMM_batch_settings; % most of the settings are here (also IPopt options)
launch_car_builder

scale.x = car.data.X_scale;
scale.u = car.data.U_scale;
scale.z = car.data.Z_scale;   


clear IPOPT_opt
IPOPT_opt = struct( ...
    'ipopt', ...
    struct( ...
    'max_iter', 5000, ... 
    'linear_solver', 'ma57', ...
    'print_level', 5, ... % 
    'ma57_pre_alloc', 10, ...
    'ma57_small_pivot_flag', 1, ...
    'ma57_pivtol', 1e-5, ...
    'mu_init', 1e-3, ...
    'warm_start_init_point', 'yes', ...
    'tol', 1e-6), ...
    'print_time', true);

%    'nlp_scaling_method', 'none', ...
    %'mu_strategy', 'adaptive', ...

Nproblems = 1; 
% Nsteps_0 = Nsteps/lap;
% Nsteps = Nsteps + 1;

% if lap > 1
%     alpha_vec_lap = linspace(0, alfa_end, (Nsteps_0 + 1));
%     alpha_vec = repmat(alpha_vec_lap(1:end-1), 1, lap-1);
%     alpha_vec = [alpha_vec, alpha_vec_lap];
% else
%     alpha_vec = linspace(0, alfa_end, Nsteps);
% end
    
%alpha_vec = linspace(0, alfa_end, Nsteps+1);

ID_instance = 1;
if init_guess 
    [init_subrange] = split_init(alpha_vec, guess, nx, 'alpha_numeric', alpha_numeric);
else
    init_subrange = [];
end

id.t{1} = [];
id.h{1} = [];
ID.t{1} = [];
ID.h{1} = [];
ID.H{1} = [];
ID.T{1} = [];
id.H{1} = [];
id.T{1} = [];
ID.start = 1;
ID.end = length(alpha_vec);
%% problem creation and solution
tic
[problem, problemData] = sub_opti_map(alpha_vec,...                                        
                                        pista,...
                                        o, id, ID,...
                                        ID_instance,...
                                        init_subrange,...
                                        alpha_vec,...
                                        IPOPT_opt,...
                                        car);
time_build = toc; % may be way faster with map class

fprintf('problem built in %f seconds \n', time_build);

activation = 0;
activation_comp = 1;
activation_opt = 1;
activation_start = 0;
Z = [];
Y = [];
RHO_head = 0;
RHO_tail = 0;

tic
solutions =  problem('x0',  problemData.w0,...
    'lbx', problemData.lbw,...
    'ubx', problemData.ubw,...
    'lbg', problemData.lbg,...
    'ubg', problemData.ubg,...
    'p',  [Z; Y; activation; RHO_head; RHO_tail; activation_comp; activation_opt; activation_start]);
time_solve = toc;

X = full(solutions.x);

%% post process
post_process;

%% Save results
problem_structure.Nproblems = 1;
problem_structure.Nsteps = Nsteps;
problem_structure.alpha_end = alfa_end;
problem_structure.alpha_vec = alpha_vec;
problem_structure.alpha_subrange = alpha_vec;
problem_structure.o = o;
problem_structure.e = e;


mkdir(savepath);
save(sprintf('%s%s',savepath, '\X_sol.mat'), 'X_sol');
save(sprintf('%s%s',savepath, '\U_sol.mat'), 'U_sol');
save(sprintf('%s%s',savepath, '\Z_sol.mat'), 'Z_sol');
save(sprintf('%s%s',savepath, '\Twist_sol.mat'), 'twist');
save(sprintf('%s%s',savepath, '\opti_var_scaled.mat'), 'X');
save(sprintf('%s%s',savepath, '\problem_structure.mat'), 'problem_structure');

%% animation

%return % comment return if you want to play the animation

animation_script;

