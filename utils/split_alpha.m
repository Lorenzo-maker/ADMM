function [alpha_subrange, id_start, id_end] = split_alpha(Nsteps, Nproblems, overlap, elemOverlap, nx, alfa_end, lap)
%
%
%
Nsteps_0 = Nsteps/lap;
Nsteps = Nsteps + 1;
if elemOverlap > nx
    overlap = overlap+1; % double states overlap (da usare solo con sequenza stati-controlli-stati)
end
% init whole alpha

if lap > 1
    alpha_range_lap = linspace(0, alfa_end, (Nsteps_0 + 1));
    alpha_range = repmat(alpha_range_lap(1:end-1), 1, lap-1);
    alpha_range = [alpha_range, alpha_range_lap];
else
    alpha_range = linspace(0, alfa_end, Nsteps);
end

% compute discretization number of each subproblem sub problems 
alpha_subsize(1:Nproblems) = floor((Nsteps)./Nproblems);
remainder = mod(Nsteps, Nproblems);
alpha_subsize(1:remainder) = alpha_subsize(1:remainder) + 1;

% split alpha in sub ranges
alpha_start = 0;
alpha_subrange = cell(Nproblems, 1);

% save indices
id_end = zeros(1, Nproblems);
id_start = zeros(1, Nproblems);

for i = 1:Nproblems
    alpha_end = alpha_start + alpha_subsize(i) + overlap;
    id_end(i) = alpha_end;
    id_start(i) = alpha_start + 1;
    % if last problem we do not have overlapping with next
    if i == Nproblems
        alpha_end = alpha_end - overlap;
    end
    alpha_subrange{i} = alpha_range(alpha_start+1: alpha_end);
    alpha_start = alpha_end - overlap;
end

end