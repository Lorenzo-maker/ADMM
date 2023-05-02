function write_consensus(ID_instance, Z, Y, RHO)
%
% function to write the consensus, ADMM multipliers, and penalty param for
% each ADMM subproblem
%

save(sprintf('Temp\\SubInstance_%i\\Z_%i.mat', ID_instance', ID_instance), 'Z');
save(sprintf('Temp\\SubInstance_%i\\Y_%i.mat', ID_instance', ID_instance), 'Y');
save(sprintf('Temp\\SubInstance_%i\\RHO_%i.mat', ID_instance', ID_instance), 'RHO');

end