% function [index_consensus, index_end] = alphaSubToIndex(alpha_subrange, overlap, Nproblems)
% index_consensus = nan(Nproblems-1, 1);
% index_end = nan(Nproblems-1, 1);
% index_consensus(1) = length(alpha_subrange{1})-floor(overlap/2)-1;
% index_end(1) = length(alpha_subrange{1});
% cumulative_sum  = length(alpha_subrange{1})- overlap;
% for i = 2:Nproblems-1
%     index_consensus(i) = cumulative_sum + length(alpha_subrange{i})-floor(overlap/2)-1;
%     index_end(i) = cumulative_sum + length(alpha_subrange{i});
%     cumulative_sum = cumulative_sum + length(alpha_subrange{i}) - overlap;
%     
% end
% index_end(end) = [];

function [index_consensus, index_end, index_consensus_local] = alphaSubToIndex(alpha_subrange, overlap, Nproblems, alpha_end, Nsteps)
index_consensus = nan(Nproblems-1, 1);
index_consensus_global = nan(Nproblems-1, 1);
index_end = nan(Nproblems-1, 1);
index_consensus(1) = (length(alpha_subrange{1})-1)-floor(overlap/2);
alfa_consesus(1) = alpha_subrange{1}(index_consensus(1));
for i = 2:Nproblems-1
    index_consensus(i) = (length(alpha_subrange{i})-1)-floor(overlap/2);
    alfa_consesus(i) = alpha_subrange{i}(index_consensus(i));
end
index_consensus_local = index_consensus;
alfa_grid = linspace(0,alpha_end,Nsteps+1);
alfa_consesus = alfa_consesus';
for i = 1:Nproblems-1
    [~,index_consensus_global(i)] = min(abs(alfa_grid-alfa_consesus(i)));
end
index_consensus = index_consensus_global;
index_end(end) = [];


