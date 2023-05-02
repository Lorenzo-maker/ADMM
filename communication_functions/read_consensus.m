function [Z, Y, RHO_head, RHO_tail] = read_consensus(ID_instance)
%
% function to read the ID subproblem consensus data updated by the ADMM
%


load(sprintf('SubInstance_%i\\Z_%i.mat', ID_instance, ID_instance), 'Z');
load(sprintf('SubInstance_%i\\Y_%i.mat', ID_instance, ID_instance), 'Y');
load(sprintf('SubInstance_%i\\RHO_%i.mat', ID_instance, ID_instance), 'RHO');

RHO_head = RHO(1);
RHO_tail = RHO(2);

end