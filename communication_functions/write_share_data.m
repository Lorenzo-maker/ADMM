function write_share_data(ID_instance, share_data)
%
% function to write the consensus, ADMM multipliers, and penalty param for
% each ADMM subproblem
%

save(sprintf('Temp\\SubInstance_%i\\share_data_%i.mat', ID_instance', ID_instance), 'share_data');


end