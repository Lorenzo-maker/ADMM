function writeSolutions(solutions, solverStats, ID_instance)

save(sprintf('Temp\\SubInstance_%i\\solutions_instance_%i.mat', ID_instance, ID_instance), 'solutions');
save(sprintf('Temp\\SubInstance_%i\\stats_instance_%i.mat', ID_instance, ID_instance), 'solverStats');

fprintf('\n');
fprintf('solutions written to file, wating for ADMM iteration...\n');
fprintf('\n');
end