function [alpha_subrange, id_start, id_end, ID, id] = split_alpha(Nsteps, Nproblems, e, o, elemOverlap, nx, alfa_end, lap)
%
%
%
Nsteps_0 = Nsteps/lap;
Nsteps = Nsteps + 1;

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
    alpha_end = alpha_start + alpha_subsize(i);
    id_end(i) = alpha_end;
    id_start(i) = max(alpha_start,1);
    % if last problem we do not have overlapping with next
    alpha_subrange{i} = alpha_range(max(alpha_start,1): alpha_end);
    alpha_start = alpha_end;
end

ID.h = cell(Nproblems,1);
ID.t = cell(Nproblems,1);
for i = 1:Nproblems
    if i == 1
        ID.h{i} = [];
        ID.t{i} = id_end(i) - floor(o/2):1:id_end(i) + ceil(o/2);
    elseif i == Nproblems
        ID.h{i} = id_start(i) - floor(o/2):1: id_start(i) + ceil(o/2);
        ID.t{i} = [];
    else
        ID.h{i} = id_start(i) - floor(o/2):1: id_start(i) + ceil(o/2);
        ID.t{i} = id_end(i) - floor(o/2):1:id_end(i) + ceil(o/2);
    end
end

ID.H = cell(Nproblems,1);
ID.T = cell(Nproblems,1);
for i = 1:Nproblems
    if i == 1
        ID.H{i} = [];
        ID.T{i} = id_end(i) + 1:1:id_end(i) + e;
    elseif i == Nproblems
        ID.H{i} = id_start(i) - e:1: id_start(i)-1;
        ID.T{i} = [];
    else
        ID.H{i} = id_start(i) - e:1: id_start(i)-1;
        ID.T{i} = id_end(i) + 1:1:id_end(i) + e;
    end
end


if e > 0 && e >= o/2 % case extended tail/head grater than consenus area
   for i = 1:Nproblems
       if i == 1
           alpha_subrange{i} = [alpha_subrange{i}, alpha_range(id_end(i) + 1: id_end(i) + e)];
       elseif i == Nproblems
           alpha_subrange{i} = [alpha_range(id_start(i) - e:id_start(i)-1), alpha_subrange{i}];
       else
           alpha_subrange{i} = [alpha_range(id_start(i) - e:id_start(i)-1), alpha_subrange{i}, alpha_range(id_end(i) + 1: id_end(i) + e)];
       end
   end
%elseif e == 0 && o > 0  % case 2: extended tail doesn't exist
else
   for i = 1:Nproblems
       if i == 1
           alpha_subrange{i} = [alpha_subrange{i}, alpha_range(id_end(i) + 1: id_end(i) + ceil(o/2))];
       elseif i == Nproblems
           alpha_subrange{i} = [alpha_range(id_start(i) - floor(o/2):id_start(i)-1), alpha_subrange{i}];
       else
           alpha_subrange{i} = [alpha_range(id_start(i) - floor(o/2):id_start(i)-1), alpha_subrange{i},alpha_range(id_end(i) + 1: id_end(i) + ceil(o/2))];
       end
   end
end

%%%% Compute index of consensus variables respect to local alpha_subrange
%%%% remember consensus are distrubuted with order right-left %%%%%%
id.h = cell(Nproblems,1);
id.t = cell(Nproblems,1);
for i = 1:Nproblems
    if i == 1
        id.h{i} = [];
        if e >= o/2 && e > 0
            id.t{i} = length(alpha_subrange{i}) - e - floor(o/2):1:length(alpha_subrange{i}) - e + ceil(o/2);
        else
            id.t{i} = length(alpha_subrange{i}) - o:1:length(alpha_subrange{i});
        end
    elseif i == Nproblems
        if e >= o/2 && e > 0
            id.h{i} = e + 1 - floor(o/2):1: e + 1 + ceil(o/2);
        else
            id.h{i} = 1:1:o+1;
        end
        id.t{i} = [];
    else
        if e >= o/2 && e > 0
            id.h{i} = e + 1 - floor(o/2):1: e + 1 + ceil(o/2);
            id.t{i} = length(alpha_subrange{i}) - e - floor(o/2):1:length(alpha_subrange{i}) - e + ceil(o/2);
        else %(caso e = 0 oppure e minore di estenzione dovuta ad o)
            id.h{i} = 1:1:o+1;
            id.t{i} = length(alpha_subrange{i}) - o:1:length(alpha_subrange{i});
        end        
    end
end

%%% index of the grid points of the extended tail/head w.r.t alpha_subrange
id.H = cell(Nproblems,1);
id.T = cell(Nproblems,1);
for i = 1:Nproblems
    if i == 1
        id.H{i} = [];
        id.T{i} = length(alpha_subrange{i}) - e + 1:1:length(alpha_subrange{i});
    elseif i == Nproblems
        id.H{i} = 1 :1: e;
        id.T{i} = [];
    else
        id.H{i} = 1 :1: e;
        id.T{i} = length(alpha_subrange{i}) - e + 1:1:length(alpha_subrange{i});
    end
end

figure(1)
ax = gca;
hold(ax,'on')
plot(1:1:Nsteps, alpha_range, 'Color', 'black', 'parent', ax)
plot(horzcat(ID.h{:}),alpha_range(horzcat(ID.h{:})),'o','Color','m','markersize',16, 'parent', ax)
plot(horzcat(ID.t{:}),alpha_range(horzcat(ID.t{:})),'*','Color','c','markersize',16, 'parent', ax)
plot(horzcat(ID.H{:}),alpha_range(horzcat(ID.H{:})),'o','Color','g', 'parent', ax)
plot(horzcat(ID.T{:}),alpha_range(horzcat(ID.T{:})),'o','Color','r', 'parent', ax)


end