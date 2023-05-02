function [solutions, solverStats] = readSolutions(ID_instance)

load(sprintf('Temp\\SubInstance_%i\\solutions_instance_%i.mat', ID_instance, ID_instance), 'solutions');
load(sprintf('Temp\\SubInstance_%i\\stats_instance_%i.mat', ID_instance, ID_instance), 'solverStats');

end