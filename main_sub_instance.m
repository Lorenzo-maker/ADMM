%% Import package and add folders

import casadi.*
% all the paths added here should be enough to run the subproblem
cd '../'
cd '../'
fprintf('\n\n %s\n\n\n', pwd)
addpath(genpath('Casadi'))
addpath(genpath('Classes'))
addpath(genpath('Track'))
addpath(sprintf('Temp\\SubInstance_%i', ID_instance));
addpath('communication_functions');
addpath(genpath('ADMM_functions'));
addpath(genpath('Data'));
addpath(genpath('scripts'));
addpath(genpath('utils'));
addpath(genpath('Extra'));

%% Create subproblems
 % load global settings
ADMM_batch_settings;

%%%%%%%%%%%%%%%%%%%%%%%% element overlap calculation %%%%%%%%%%%%%%%%%%%%%%
if overlap >= 1 && overlap_control == true
    elemOverlap = (overlap*(nx + nu + nz)+nx);
else
    elemOverlap = nx;
end

% split alpha
[alpha_subrange] = split_alpha(Nsteps, Nproblems, overlap, elemOverlap, nx, alfa_end, lap);
% what is init_subrange used for ?
alpha_subrange = alpha_subrange{ID_instance};
[init_subrange] = split_init(alpha_subrange, guess, nx+1, nu);

convergence = 0;
problem = [];
problemData = [];
sol_struct = [];

tic

waiting = true;
waiting_filename = sprintf('Temp\\%s_%i\\%s_%i.txt', 'SubInstance', ID_instance, 'waiting_file', ID_instance);

if ID_instance == 1
    overlap_tail = overlap;
    overlap_head = 0;
elseif ID_instance == Nproblems
    overlap_tail = 0;
    overlap_head = overlap;
else
    overlap_tail = overlap;
    overlap_head = overlap;
end

% pista = track_fun(1, '3D', 1, 'end', 1500, 'alfa_limited', true, 'alfa_lim', [max(0, alpha_subrange(1)-10*dalfa), min(alfa_end, alpha_subrange(end)+10*dalfa)]);
load(sprintf('Temp//SubInstance_%i//pista.mat',ID_instance))
[problem, problemData] = sub_opti_map(alpha_subrange,...
                                        pista,...
                                        overlap_tail, overlap_head,...
                                        ID_instance, direct_c,...
                                        init_subrange,...
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

%     elseif ADMM_iteration == 3
%         activation = 0;
%         activation_comp = 10^2;
%         activation_opt = 1;
%         activation_start = 0;
% 
%         solutions =  problem('x0',  problemData.w0,...
%             'lbx', problemData.lbw,...
%             'ubx', problemData.ubw,...
%             'lbg', problemData.lbg,...
%             'ubg', problemData.ubg,...
%             'p',  [Z; Y; activation; RHO_head; RHO_tail; activation_comp; activation_opt; activation_start],...
%             'lam_g0', full(sol_struct.lam_g),...
%             'lam_x0', full(sol_struct.lam_x));        
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