%% ADMM main instance script
clear all
clc
close all
working_dir = pwd;

savepath = 'Results\1500_4_81'; % for results path

%% Path

addpath('communication_functions');
addpath(genpath('Casadi'));
addpath(genpath('Classes'));
addpath(genpath('Track'));
addpath(genpath('ADMM_functions'));
addpath(genpath('Data'));
addpath(genpath('scripts'));
addpath(genpath('utils'));
addpath(genpath('Extra'));

import casadi.*

%% Create track
% pista = track_fun(1, '3D', 1, 'end', 700,'degree',4);
warning off
load('Data\pista_original_2.mat');
warning on
pista.residual;

%pista = track2casadiSpline(pista);

disp(['fitting error is ', num2str(full(pista.tot_error))])

save('Data\pista.mat', 'pista');

%% load ADMM settings
ADMM_batch_settings; % most of the settings are here (also IPopt options)

if o > 0 
    elemOverlap = (o*(nx + nu + nz)+nx);
else
    elemOverlap = nx;
end

J0 = zeros(1, Nproblems);
[alpha_subrange] = split_alpha(Nsteps, Nproblems, e, o, elemOverlap, nx, alfa_end, lap);

% create problems
problem = cell(Nproblems, 1);
solver_info = cell(Nproblems, 1);
problemData = cell(Nproblems, 1);
sol_struct = cell(Nproblems, 1);
RHO_head = 1.*ones(Nproblems, 1); % penalty head
RHO_head_0 = 1.*ones(Nproblems, 1); % penalty head 0
RHO_tail = 1.*ones(Nproblems, 1); % penalty tail
RHO_tail_0 = 1.*ones(Nproblems, 1); % penalty tail 0
previousConsErr = 0;

N_OCP = 0;
for i = 1:Nproblems

    alpha_subrange{i} = round(alpha_subrange{i}, 15); % elemOverlap need to be > nx
    
end

% definite another time effective element overlap 
elemOverlap = (1*(nx + nu + nz) + nx);

%% plots init

ADMM_plot_init;
figADMM = gcf;

%% ADMM_algorithm
% variable initialization
solutions = cell(Nproblems, 1);
convergence = false;
X = cell(1, Nproblems); % init solution arrays for each sub problem
Xc = cell(1, Nproblems); % init solution arrays for each colloc
Z = cell(1, Nproblems); % consensus for sub - problem
Z_acc = cell(1, Nproblems); % consensus acc for sub - problem
Z_acc_all = cell(1, Nproblems); % consensus acc for all problem 
Y = cell(1, Nproblems); % moltiplicatori sub - problem
Y_acc = cell(1, Nproblems); % moltiplicatori acc sub - problem
Y_acc_all = cell(1, Nproblems); % moltiplicatori acc for all problem 
gamma = cell(1, Nproblems); % gamma for all problem
error = cell(1, Nproblems);
Znext = cell(1, Nproblems);
Ynext = cell(1, Nproblems);
consensusChange_norm = cell(1, Nproblems);
consensusChange = cell(1, Nproblems);
consErr_norm = cell(1, Nproblems);
consErr = cell(1, Nproblems);
conv = cell(1, Nproblems);
r_head = nan(1, Nproblems);
r_tail = nan(1, Nproblems);
s_head = nan(1, Nproblems);
s_tail = nan(1, Nproblems);
epsilon_error = cell(1, Nproblems);
epsilon_z = cell(1, Nproblems);
waiting_filename = cell(1, Nproblems);
        
% create consensus arrays multipliers and tollerance array
for i = 1:Nproblems
    if i == 1 % no head 
        epsilon_error{i} = error_factor*epsilon_array;
        epsilon_z{i} = delta_z_factor*epsilon_array;
        numConsElem = elemOverlap;
    elseif i == Nproblems  % no tail
        epsilon_error{i} = error_factor*epsilon_array;
        epsilon_z{i} = delta_z_factor*epsilon_array;
        numConsElem = elemOverlap;
    else % head and tail
        epsilon_error{i} = error_factor*[epsilon_array;epsilon_array];
        epsilon_z{i} = delta_z_factor*[epsilon_array;epsilon_array];
        numConsElem = 2*elemOverlap; 
    end
    Z{i} = zeros(numConsElem, 1);      % consensus parameters
    X{i} = 0;                    % preallocate cells also for solutions (will be overwritten after iteration 0)
    Y{i} = full(Z{i}.*0);              % consensus multipliers
    conv{i} = 0;
end
Z00 = Z;
Y00 = Y;
% Initialize other quantities
Zprevious = Z;
Yprevious = Y;
ADMM_iteration = 0;
check_gamma = [];

%% SUB-INSTANCES INITIALIZATION
%write initial data to files, convergence, Z, Y, RHO
parfor i = 1:Nproblems

    mkdir(sprintf('%s\\Temp\\SubInstance_%i', working_dir, i));
    copyfile('Data\pista.mat',sprintf('Temp\\SubInstance_%i\\pista.mat', i));
    waiting_filename{i} = create_waiting_function(i);
    write_convergence(i, convergence);
    write_consensus(i, Z{i}, Y{i}, [RHO_head(i); RHO_tail(i)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% HERE WE LAUNCH INSTANCES IN BATCH %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pause(i*1);
    launch_sub_instance(i, working_dir);
end


% recap: waiting status = 0 -> problems are working
%        waiting status = 1 -> problems are completed, ADMM can gather
%                              solutions
waiting_status = 0;
waiting = nan(Nproblems, 1);
while ~waiting_status % loop until all subrpoblems are are ready
    
    for i = 1: Nproblems
        waiting(i) = readWaiting(waiting_filename{i});
    end
    
    waiting_status = all(waiting);
    pause(5);
    % when all the problems turned to 1 (casadi problems built), ADMM can
    % start the sub problems iterations
end

%%%% make the problems finally start
for i = 1:Nproblems
   change_waiting_status(i); % turning waiting to -> 0
end
pause(3);

%% ADMM LOOP
tic
while ~convergence % test di convergenza del consenso
    bool_out = [];
    %%%%%%% consensus activation parameter is handled directly inside each 
    %%%%%%% subproblem
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% need to wait the subproblems to complete
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% by reading the wating status file
    % recap: waiting status = 0 -> problems are working
    %        waiting status = 1 -> problems are completed, ADMM can gather
    %                              solutions
    waiting_status = 0;
    waiting = nan(Nproblems, 1);
    while ~waiting_status % loop until all subrpoblems are are ready
        
        for i = 1: Nproblems
            wait = [];
            counter = 0;
            while isempty(wait)
                wait = readWaiting(waiting_filename{i});
                counter = counter+1;
                if counter > 1
                    disp('error reading waiting file... trying again')
                    pause(1);
                end
            end
            waiting(i) = wait;
        end
        
        waiting_status = all(waiting);
        pause(3);
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD SOLUTIONS
    parfor kk = 1:Nproblems
        
         [solutions{kk}, solver_info{kk}] = readSolutions(kk);
        
        X{kk} = full(solutions{kk}.x);
        if direct_c
            %reshape solution to delete collocation points
            %%%%% assumption d is equal for every subproblem for now...%%%%
            X{kk} = [];
            X_origin{kk} = full(solutions{kk}.x);
            X{kk} = reshape(X_origin{kk}(nx+1:end), [], length(alpha_subrange{kk})-1);
            Xc{kk} = X{kk}(nu+nz+1:nu+nz+d*nx,:);                       % (nu+nz+1:nu+nz+problemData{kk}.d*nx,:)
            X{kk} = [X{kk}(1:nu+nz,:);X{kk}(nu+nz+d*nx+1:end,:)];
            X{kk} = X{kk}(:);
            X{kk} = [X_origin{kk}(1:nx);X{kk}];
        end
        
        
        sol_struct{kk} = solutions{kk};
        
        if ADMM_iteration == 2 %1
            J0(kk) = full(solutions{kk}.f);
        end        
        
    end
    
    solver_info_admm{ADMM_iteration+1} = solver_info;
    
    % aggiornamento dei parametri
    for i = 1:Nproblems
        
        % aggiorno il consenso locale
        if i == 1 % no head
            
            X_tail = X{i}(end - (nx + (floor(overlap/2) + 1)*(nx + nu + nz))+1: end - (nx + floor(overlap/2 + 1)*(nx + nu + nz))+(nx+nu+nz+nx));
            X_head_next = X{i+1}(floor(overlap/2)*(nx+nu+nz)+1: floor(overlap/2)*(nx+nu+nz)+(nx+nu+nz+nx));
            Y_tail = Y{i};
            if Nproblems > 2
                Y_head_next = Y{i+1}(1: length(Y{i+1})/2);
            else
                Y_head_next = Y{i+1};
            end                
            Ztail = update_dual_var(X_tail, X_head_next, RHO_tail(i), RHO_head(i+1), Y_tail, Y_head_next, 'method', 'Augmented-Lagrangian');
            
            X_overlap = X_tail;
            error{i} = (X_overlap - Z{i});
            
            Z{i} = Ztail; %Z barra
            Z_acc{i} = update_var_acc(Z{i}, Zprevious{i} , Z_acc_all, i, ADMM_iteration, 'method', update_method,'alfa', alfa_stationary); %Z accelerato
            
            Y{i} = Y{i} + RHO_tail(i).*(X_overlap - Z{i}); % Y barra
            Y_acc{i} = update_var_acc(Y{i}, Yprevious{i} , Y_acc_all, i, ADMM_iteration, 'method', update_method,'alfa', alfa_stationary); %Y accelerato
            
            % aggiorno l'errore inteso come differenza della soluzione di
            % sovrapposizione corrente dal consenso corrente
            %error{i} = (X_overlap - Z{i});
            
            if ADMM_iteration == 2
                RHO_tail(i) = J0(i).*2./vecnorm(X_overlap - Zprevious{i}).^2/rho_scale_tail; % error{i}
                RHO_tail_0(i) = J0(i).*2./vecnorm(X_overlap - Zprevious{i}).^2/rho_scale_tail; % error{i}
            end
            
            if all(conv{i}) == false && ADMM_iteration > 1
                RHO_tail(i) = updateRHO(RHO_tail(i), X_overlap, Z{i}, Zprevious{i}(end - elemOverlap + 1 : end),  Y{i}(end - elemOverlap + 1 : end), ADMM_iteration);
            end
            
            % aggiorno i moltiplicatori locali

            r_head(i) = 0;
            r_tail(i) = vecnorm(error{i});
            s_head(i) = 0;
            s_tail(i) = vecnorm((Z{i} - Zprevious{i}));
            
        elseif i == Nproblems % no tail
            
            X_overlap = X{i}(floor(overlap/2)*(nx+nu+nz)+1: floor(overlap/2)*(nx+nu+nz)+(nx+nu+nz+nx));            
            error{i} = (X_overlap - Z{i});
            
            Z{i} = Ztail; % Z barra
            Z_acc{i} = update_var_acc(Z{i}, Zprevious{i}, Z_acc_all, i, ADMM_iteration, 'method', update_method,'alfa', alfa_stationary); %Z accelerato
            
            Y{i} = Y{i} + RHO_head(i).*(X_overlap - Z{i});% Y barra
            Y_acc{i} = update_var_acc(Y{i}, Yprevious{i} , Y_acc_all, i, ADMM_iteration, 'method', update_method,'alfa', alfa_stationary); %Y accelerato
            
            % aggiorno l'errore inteso come differenza della soluzione di
            % sovrapposizione corrente dal consenso corrente
            %error{i} = (X_overlap - Z{i});
            
            if ADMM_iteration == 2
                RHO_head(i) = J0(i).*2./vecnorm(X_overlap - Zprevious{i}).^2/rho_scale_head;%/4   error{i}
                RHO_head_0(i) = J0(i).*2./vecnorm(X_overlap - Zprevious{i}).^2/rho_scale_head;%/4  error{i}                 
            end
            
            if all(conv{i}) == false && ADMM_iteration > 1
                RHO_head(i) = updateRHO(RHO_head(i), X_overlap, Z{i}, Zprevious{i}(1:elemOverlap),  Y{i}(1:elemOverlap), ADMM_iteration);
            end
            
            % aggiorno i moltiplicatori locali

            
            r_head(i) = vecnorm(error{i});
            r_tail(i) = 0;
            s_head(i) = vecnorm((Z{i} - Zprevious{i}));
            s_tail(i) = 0;
            
        else % tail and head
            
            Zhead = Ztail; % share with the previous tail problem (Z_i-1)
            
            % update Z_i
            X_tail = X{i}(end - (nx + floor(overlap/2 + 1)*(nx + nu + nz))+1: end - (nx + floor(overlap/2 + 1)*(nx + nu + nz))+(nx+nu+nz+nx));
            X_head_next = X{i+1}(floor(overlap/2)*(nx+nu+nz)+1: floor(overlap/2)*(nx+nu+nz)+(nx+nu+nz+nx));
            Y_tail = Y{i}(length(Y{i})/2 + 1:end);
            Y_head_next = Y{i+1}(1: elemOverlap);
            Ztail = update_dual_var(X_tail, X_head_next, RHO_tail(i), RHO_head(i+1), Y_tail, Y_head_next, 'method', 'Augmented-Lagrangian');
            
            X_head = X{i}(floor(overlap/2)*(nx+nu+nz)+1: floor(overlap/2)*(nx+nu+nz)+(nx+nu+nz+nx));
            X_overlap = [X_head; X_tail];            
            error{i} = (X_overlap - Z{i});
            
            Z{i} = [Zhead; Ztail]; % Z barra
            Z_acc{i} = update_var_acc(Z{i}, Zprevious{i} , Z_acc_all, i, ADMM_iteration, 'method', update_method,'alfa', alfa_stationary); %Z accelerato
            
            Y{i} = Y{i} + [RHO_head(i).*(X_head - Zhead); RHO_tail(i).*(X_tail - Ztail)]; %Y barra
            Y_acc{i} = update_var_acc(Y{i}, Yprevious{i} , Y_acc_all, i, ADMM_iteration, 'method', update_method,'alfa', alfa_stationary); %Y accelerato
            
            
            % aggiorno l'errore inteso come differenza della soluzione di
            % sovrapposizione corrente dal consenso corrente
            %error{i} = (X_overlap - Z{i});
            
            if ADMM_iteration == 2
                RHO_tail(i) = J0(i).*1./vecnorm(X_overlap(end/2+1:end) - Zprevious{i}(end - elemOverlap + 1 : end)).^2/rho_scale_tail; % error{i}
                RHO_tail_0(i) = J0(i).*1./vecnorm(X_overlap(end/2+1:end) - Zprevious{i}(end - elemOverlap + 1 : end)).^2/rho_scale_tail; % error{i}
                RHO_head(i) = J0(i).*1./vecnorm(X_overlap(1:end/2) - Zprevious{i}(1:elemOverlap)).^2/rho_scale_head; % error{i}
                RHO_head_0(i) = J0(i).*1./vecnorm(X_overlap(1:end/2) - Zprevious{i}(1:elemOverlap)).^2/rho_scale_head; % error{i}                
            end
            
            if all(conv{i}) == false && ADMM_iteration > 1
                RHO_head(i) = updateRHO(RHO_head(i), X_head, Zhead, Zprevious{i}(1:elemOverlap), Y{i}(1:elemOverlap), ADMM_iteration);
                RHO_tail(i) = updateRHO(RHO_tail(i), X_tail, Ztail, Zprevious{i}(end - elemOverlap + 1 : end), Y{i}(end - elemOverlap + 1 : end), ADMM_iteration);
            end
            
            % aggiorno i moltiplicatori locali

            r_head(i) = vecnorm(error{i}(1:end/2));
            r_tail(i) = vecnorm(error{i}(end/2+1:end));
            s_head(i) = vecnorm((Z{i}(1:end/2) - Zprevious{i}(1:elemOverlap)));
            s_tail(i) = vecnorm((Z{i}(end/2+1:end) - Zprevious{i}(end - elemOverlap + 1 : end)));
        end
        
    end
    
    % Memorizzo Z e Y accelerati
    Z_acc_all{ADMM_iteration + 1} = Z_acc;
    Y_acc_all{ADMM_iteration + 1} = Y_acc;
    
    % Aggiornamento del gamma
    RHOx = mean([RHO_head;RHO_tail]);
    gamma{ADMM_iteration + 1} = (1/RHOx)*vecnorm(vertcat(Y_acc{:})-vertcat(Yprevious{:}))^2 + RHOx*vecnorm(vertcat(Z_acc{:})-vertcat(Zprevious{:}))^2;
    if ADMM_iteration < 1
        gamma0 = 0;
    else
        gamma0 = gamma{1};
    end
    
    % Condizione per aggiornamento Z e Y
    if gamma{ADMM_iteration + 1} < gamma0*eta^(ADMM_iteration+1)
        Znext = Z_acc;
        Ynext = Y_acc;
        check_gamma = [check_gamma; 1];
    else
        Znext = Z;
        Ynext = Y;
        check_gamma = [check_gamma; 0];
        gamma{ADMM_iteration + 1} = (1/RHOx)*vecnorm(vertcat(Ynext{:})-vertcat(Yprevious{:}))^2 + RHOx*vecnorm(vertcat(Znext{:})-vertcat(Zprevious{:}))^2;
    end

    % To avoid new Z and Y after the preliminary optimization
    if ADMM_iteration < 2
        Z = Z00;
        Y = Y00;
    elseif ADMM_iteration == 2
        Z = Znext;
        Y = Y00;
    else
        Z = Znext;
        Y = Ynext;
    end
    
    % Verifica convergenza
    for i = 1:Nproblems
        consensusChange_norm{i} = vecnorm((Znext{i} - Zprevious{i})); % norm of dual residual
        consensusChange{i} = abs(Znext{i} - Zprevious{i}); % dual residual
        consErr_norm{i} = vecnorm(error{i}); % norm of primal residual
        consErr{i} = abs(error{i}); % primal residual
        %conv_old{i} = (consErr_eu{i} < consensusTol.*sqrt(length(Znext{i}))) & (consensusChange_eu{i} < consensusChangeTol.*sqrt(length(Znext{i})));
        conv{i} = (consErr{i} < epsilon_error{i}) & (consensusChange{i} < epsilon_z{i});
        bool_out = [bool_out; conv{i}];
    end
    
    % test di convergenza
    convergence = all(bool_out); % check if all are true
    
    % output printing
    ADMM_print(ADMM_iteration, Nproblems, consErr_norm, consensusChange_norm, RHO_head, RHO_tail, conv, solver_info);
    

    
    % aggiorno il vettore dei consensi
    Zprevious = Znext;
    Yprevious = Ynext;
    previousConsErr = consErr;
    ADMM_iteration = ADMM_iteration + 1;
    Xinter{ADMM_iteration} = X;
    Zinter{ADMM_iteration} = Z;
    
    % writing data to subproblems
    parfor i = 1:Nproblems
        % maybe just 1 function to write all data togheder to a structure is enough ??
        write_convergence(i, convergence);
        write_consensus(i, Z{i}, Y{i}, [RHO_head(i); RHO_tail(i)]);
        % when data is written just update the waiting status to launch the
        % next subproblems iteration
        pause(1);
        change_waiting_status(i);
    end
    
    % plots update
    ADMM_plot_update;

end
time_while = toc;
delete(gcp);

%% post processing whole solution (plots)
post_process;

%% Save results
problem_structure.Nproblems = Nproblems;
problem_structure.Nsteps = Nsteps;
problem_structure.alpha_end = alfa_end;
problem_structure.alpha_vec = alpha_vec;
problem_structure.alpha_subrange = alpha_subrange;
problem_structure.overlap = overlap;


mkdir(savepath);
savefig(figADMM,[savepath, '\ADMM_plot'], 'compact');
close(figADMM);
save(sprintf('%s%s',savepath, '\X_sol.mat'), 'X_sol');
save(sprintf('%s%s',savepath, '\U_sol.mat'), 'U_sol');
save(sprintf('%s%s',savepath, '\Z_sol.mat'), 'Z_sol');
save(sprintf('%s%s',savepath, '\Twist_sol.mat'), 'twist');
save(sprintf('%s%s',savepath, '\opti_var_scaled.mat'), 'X');
save(sprintf('%s%s',savepath, '\problem_structure.mat'), 'problem_structure');

%% gif animation (not working)
% gif_creator_script;

return % comment return if you want to play the animation

animation_script;