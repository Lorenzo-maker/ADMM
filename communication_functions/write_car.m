function write_car(ID_instance, car)
%
% function to write the consensus, ADMM multipliers, and penalty param for
% each ADMM subproblem
%

save(sprintf('Temp\\SubInstance_%i\\car_%i.mat', ID_instance', ID_instance), 'car');


end