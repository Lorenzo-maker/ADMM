function ADMM_print(ADMM_iteration, Nproblems, consErr, consensusChange, RHO_h, RHO_t, convergence, solver_info)

charlen = length(num2str(ADMM_iteration));
iter_array = repmat('_', Nproblems, charlen);
iter_array(1,:) = num2str(ADMM_iteration);
fprintf('\n******************************************************************************************\n')
fprintf('it%i     ||r||       ||d_Z||           rho_h           rho_t         converg         sub_conv         t_nlp         iter\n', ADMM_iteration);
error_array = cat(2, [consErr{:}]);
change_array = cat(2, [consensusChange{:}]);
for i = 1:Nproblems
    fprintf('%i   %e  %e    %e    %e         %i               %i           %f        %i\n', i, error_array(i), change_array(i), RHO_h(i), RHO_t(i), all(convergence{i}), solver_info{i}.success, solver_info{i}.t_proc_total, solver_info{i}.iter_count); % iter_array(i,:)
end
fprintf('*********************************************************************************************\n')
end