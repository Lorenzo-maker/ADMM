%% Import package and add folders

import casadi.*
% all the paths added here should be enough to run the subproblem
cd '../'
cd '../'
fprintf('\n\n %s\n\n\n', pwd)
addpath(genpath('../../../Casadi'))
addpath(genpath('Classes'))
addpath(genpath('Track'))
addpath(sprintf('Temp\\SubInstance_%i', ID_instance));
addpath('communication_functions');
addpath(genpath('ADMM_functions'));
addpath(genpath('Data'));
addpath(genpath('scripts'));
addpath(genpath('utils'));

%% Create subproblems
 % load global settings
load(sprintf('Temp//SubInstance_%i//pista.mat',ID_instance))
ADMM_batch_settings;
if split_manual
    load(sprintf('Temp//SubInstance_%i//manual_index.mat',ID_instance))
else
    manual_index = [];
end
%%%%%%%%%%%%%%%%%%%%%%%% element overlap calculation %%%%%%%%%%%%%%%%%%%%%%

elemOverlap = (o*(nx + nu + nz)+nx);

%%%%% Load alpha_subrange
[alpha_subrange, ~, ~, ID, id] = split_alpha(Nsteps, Nproblems, e, o, alpha_vec, lap, manual_index);
alpha_subrange = alpha_subrange{ID_instance};
%%%%% Use initial guess from simulation
if init_guess 
    [init_subrange] = split_init(alpha_subrange, guess, nx, 'alpha_numeric', alpha_numeric);
else
    init_subrange = [];
end

convergence = 0;
problem = [];
problemData = [];
sol_struct = [];

tic

waiting = true;
waiting_filename = sprintf('Temp\\%s_%i\\%s_%i.txt', 'SubInstance', ID_instance, 'waiting_file', ID_instance);

% if ID_instance == 1
%     overlap_tail = e;
%     overlap_head = 0;
% elseif ID_instance == Nproblems
%     overlap_tail = 0;
%     overlap_head = e;
% else
%     overlap_tail = e;
%     overlap_head = e;
% end
[problem, problemData] = sub_opti_map(alpha_subrange,...                                        
                                        pista,...
                                        o, id, ID,...
                                        ID_instance, d,...
                                        init_subrange,...
                                        alpha_vec,...
                                        IPOPT_opt);
% save(sprintf('Temp//SubInstance_%i//longData.mat',ID_instance), 'problemData');
fprintf('\n problem %i built \n', ID_instance);
change_waiting_status(ID_instance); % change here to 1 to let ADMM know the problems are built

toc


ADMM_iteration = 0;
while ~convergence
    

    while waiting == true % wait for ADMM update
        wait = [];
        counter = 0;
        while isempty(wait)
            wait = readWaiting(waiting_filename);
            counter = counter+1;
            if counter > 1
                disp('error reading waiting file... trying again')
                pause(1);
            end
        end
        waiting = wait;
        pause(3); % wait 5 seconds between each file reading
    end
    fprintf('problem_number %i computing_iteration %i \n', ID_instance, ADMM_iteration)
    convergence = read_convergence(ID_instance);
    if convergence
%         return % return
        exit()  % directly close the program if convergence is achieved 
    end
    
    [Z, Y, RHO_head, RHO_tail] = read_consensus(ID_instance);
    
%     save(sprintf('Z_%i_%i.mat', ID_instance, ADMM_iteration), 'Z');
%     save(sprintf('Y_%i_%i.mat', ID_instance, ADMM_iteration), 'Y');
    %%% waiting = 0 -> start computations...
    if homotopy
        if ADMM_iteration == 0
            fprintf('computing iteration 0')
            activation = 0;
            activation_comp = 10^-2;%1;
            activation_opt = 0;
            activation_start = 1;

            solutions =  problem('x0',  problemData.w0,...
                'lbx', problemData.lbw,...
                'ubx', problemData.ubw,...
                'lbg', problemData.lbg,...
                'ubg', problemData.ubg,...
                'p',  [Z; Y; activation; RHO_head; RHO_tail; activation_comp; activation_opt; activation_start]);

        elseif ADMM_iteration == 1
            activation = 0;
            activation_comp = 10^-1;%10^2;
            activation_opt = 0.01;%0.01;
            activation_start = 0.99;%0.8;

            solutions =  problem('x0',  problemData.w0,...
                'lbx', problemData.lbw,...
                'ubx', problemData.ubw,...
                'lbg', problemData.lbg,...
                'ubg', problemData.ubg,...
                'p',  [Z; Y; activation; RHO_head; RHO_tail; activation_comp; activation_opt; activation_start],...
                'lam_g0', full(sol_struct.lam_g),...
                'lam_x0', full(sol_struct.lam_x));

         %%%...if necessary add some other homotopy iteration...%%%   
        elseif ADMM_iteration == 2
            activation = 0;
            activation_comp = 1;%10^2;
            activation_opt = 1;
            activation_start = 0;

            solutions =  problem('x0',  problemData.w0,...
                'lbx', problemData.lbw,...
                'ubx', problemData.ubw,...
                'lbg', problemData.lbg,...
                'ubg', problemData.ubg,...
                'p',  [Z; Y; activation; RHO_head; RHO_tail; activation_comp; activation_opt; activation_start],...
                'lam_g0', full(sol_struct.lam_g),...
                'lam_x0', full(sol_struct.lam_x));       
        else
            activation = 1;
            activation_comp = 1;%10^2;
            activation_opt = 1;
            activation_start = 0;

            solutions =  problem('x0',  problemData.w0,...
                'lbx', problemData.lbw,...
                'ubx', problemData.ubw,...
                'lbg', problemData.lbg,...
                'ubg', problemData.ubg,...
                'p',  [Z; Y; activation; RHO_head; RHO_tail; activation_comp; activation_opt; activation_start],...
                'lam_g0', full(sol_struct.lam_g),...
                'lam_x0', full(sol_struct.lam_x));
        end
    else
        if ADMM_iteration == 0
            activation = 1;
            activation_comp = 1;%10^2;
            activation_opt = 1;
            activation_start = 0;

            solutions =  problem('x0',  problemData.w0,...
                'lbx', problemData.lbw,...
                'ubx', problemData.ubw,...
                'lbg', problemData.lbg,...
                'ubg', problemData.ubg,...
                'p',  [Z; Y; activation; RHO_head; RHO_tail; activation_comp; activation_opt; activation_start]);
        else
            activation = 1;
            activation_comp = 1;%10^2;
            activation_opt = 1;
            activation_start = 0;

            solutions =  problem('x0',  problemData.w0,...
                'lbx', problemData.lbw,...
                'ubx', problemData.ubw,...
                'lbg', problemData.lbg,...
                'ubg', problemData.ubg,...
                'p',  [Z; Y; activation; RHO_head; RHO_tail; activation_comp; activation_opt; activation_start],...
                'lam_g0', full(sol_struct.lam_g),...
                'lam_x0', full(sol_struct.lam_x));
        end
    end
        
    %%%
    
    % save files ...
    
    %%%
    sol_struct = solutions;
    problemData.w0 = full(solutions.x);
    solutions.f = full(solutions.f);
    solutions.g = full(solutions.g);
    solutions.lam_g = full(solutions.lam_g);
    solutions.lam_p = full(solutions.lam_p);
    solutions.lam_x = full(solutions.lam_x);
    solutions.x = full(solutions.x);
    fprintf('problem %i finished...\n', ID_instance);
    writeSolutions(solutions, problem.stats(), ID_instance);
    pause(2); % make a pause to be sure all data had the time to be properly written and set in the directories
    change_waiting_status(ID_instance); % switching waiting to -> 0 (ADMM can start)
    pause(0.5)
    waiting = readWaiting(waiting_filename);

    ADMM_iteration = ADMM_iteration+1;
end

exit(); % exit subinstance command window (needed when parallelizing with many windows, can be commented for debugging)