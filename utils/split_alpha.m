function [alpha_subrange, id_start, id_end, ID, id] = split_alpha(Nsteps, Nproblems, e, o, alpha_range, lap, manual_index)
%
%
%
Nsteps_0 = Nsteps/lap;
Nsteps = Nsteps + 1;
% if lap > 1
%     alpha_range_lap = linspace(0, alfa_end, (Nsteps_0 + 1));
%     alpha_range = repmat(alpha_range_lap(1:end-1), 1, lap-1);
%     alpha_range = [alpha_range, alpha_range_lap];
% else
%     alpha_range = linspace(0, alfa_end, Nsteps);
% end

% compute discretization number of each subproblem sub problems 
if isempty(manual_index)
    alpha_subsize(1:Nproblems) = floor((Nsteps)./Nproblems);
    remainder = mod(Nsteps, Nproblems);
    alpha_subsize(1:remainder) = alpha_subsize(1:remainder) + 1;
else
    if lap == 1
        manual_index = [0, manual_index, length(alpha_range)];
        alpha_subsize = diff(manual_index);
    else
        manual_index = [0, manual_index, (length(alpha_range)-1)/lap];
        alpha_subsize = repmat(diff(manual_index), 1, lap);
        alpha_subsize(end) = alpha_subsize(end) + (length(alpha_range) - sum(alpha_subsize));
    end
end

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

% Store index of consensus mesh interval relative to the whole alphagrid
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

% Store index of extended H/T mesh interval relative to the whole alphagrid
ID.H = cell(Nproblems,1);
ID.T = cell(Nproblems,1);
for i = 1:Nproblems
    if i == 1
        ID.H{i} = [];
        ID.T{i} = id_end(i) - floor(o/2) + 1 :1:min(id_end(i) - floor(o/2) + e, Nsteps);
        %ID.T{i} = id_end(i) - floor(o/2) + 1 :1:id_end(i) - floor(o/2) + e;
    elseif i == Nproblems
        ID.H{i} = max(1, id_start(i) + ceil(o/2) - e):1: id_start(i) + ceil(o/2) - 1;
        %ID.H{i} = id_start(i) + ceil(o/2) - e:1: id_start(i) + ceil(o/2) - 1;
        ID.T{i} = [];
    else
        ID.H{i} = max(1, id_start(i) + ceil(o/2) - e):1: id_start(i) + ceil(o/2) - 1;
        %ID.H{i} = id_start(i) + ceil(o/2) - e:1: id_start(i) + ceil(o/2) - 1;
        ID.T{i} = id_end(i) - floor(o/2) + 1 :1:min(id_end(i) - floor(o/2) + e, Nsteps);
        %ID.T{i} = id_end(i) - floor(o/2) + 1 :1:id_end(i) - floor(o/2) + e;
    end
end

for i = 1:Nproblems
   if i == 1
       %alpha_subrange{i} = [alpha_subrange{i}, alpha_range(id_end(i) + 1: id_end(i) + ceil(o/2) + max(0, e-o))];
       alpha_subrange{i} = [alpha_subrange{i}, alpha_range(id_end(i) + 1: min(id_end(i) + ceil(o/2) + max(0, e-o), Nsteps))];
   elseif i == Nproblems
       %alpha_subrange{i} = [alpha_range(id_start(i) - floor(o/2) - max(0, e-o):id_start(i)-1), alpha_subrange{i}];
       alpha_subrange{i} = [alpha_range(max(1, id_start(i) - floor(o/2) - max(0, e-o)):id_start(i)-1), alpha_subrange{i}];
   else
       %alpha_subrange{i} = [alpha_range(id_start(i) - floor(o/2) - max(0, e-o):id_start(i)-1), alpha_subrange{i},alpha_range(id_end(i) + 1: id_end(i) + ceil(o/2) + max(0, e-o))];
       alpha_subrange{i} = [alpha_range(max(1,id_start(i) - floor(o/2) - max(0, e-o)):id_start(i)-1), alpha_subrange{i}, alpha_range(id_end(i) + 1: min(id_end(i) + ceil(o/2) + max(0, e-o),Nsteps))];
   end
end




% if e > 0 && e > o % case extended tail/head grater than consenus area
%    for i = 1:Nproblems
%        if i == 1
%            alpha_subrange{i} = [alpha_subrange{i}, alpha_range(id_end(i) + 1: id_end(i) + et)];
%        elseif i == Nproblems
%            alpha_subrange{i} = [alpha_range(id_start(i) - eh:id_start(i)-1), alpha_subrange{i}];
%        else
%            alpha_subrange{i} = [alpha_range(id_start(i) - eh:id_start(i)-1), alpha_subrange{i}, alpha_range(id_end(i) + 1: id_end(i) + et)];
%        end
%    end
% %elseif e == 0 && o > 0  % case 2: extended tail doesn't exist
% else
%    for i = 1:Nproblems
%        if i == 1
%            alpha_subrange{i} = [alpha_subrange{i}, alpha_range(id_end(i) + 1: id_end(i) + ceil(o/2))];
%        elseif i == Nproblems
%            alpha_subrange{i} = [alpha_range(id_start(i) - floor(o/2):id_start(i)-1), alpha_subrange{i}];
%        else
%            alpha_subrange{i} = [alpha_range(id_start(i) - floor(o/2):id_start(i)-1), alpha_subrange{i},alpha_range(id_end(i) + 1: id_end(i) + ceil(o/2))];
%        end
%    end
% end

%%%% Compute index of consensus variables respect to local alpha_subrange
%%%% remember consensus are distrubuted with order right-left %%%%%%
id.h = cell(Nproblems,1);
id.t = cell(Nproblems,1);
for i = 1:Nproblems
    if i == 1
        id.h{i} = [];
        id.t{i} = length(alpha_subrange{i}) - max(e,o):1:length(alpha_subrange{i}) - max(e,o) + o;
    elseif i == Nproblems
        id.h{i} = 1 + max(e, o) - o:1: 1 + max(e, o);
        id.t{i} = [];
    else
        id.t{i} = length(alpha_subrange{i}) - max(e,o):1:length(alpha_subrange{i}) - max(e,o) + o;
        id.h{i} = 1 + max(e, o) - o:1: 1 + max(e, o);
    end
end

%%% index of the grid points of the extended tail/head w.r.t alpha_subrange
id.H = cell(Nproblems,1);
id.T = cell(Nproblems,1);
for i = 1:Nproblems
    if i == 1
        id.H{i} = [];
        id.T{i} = length(alpha_subrange{i}) - max(e, o) + 1:1:length(alpha_subrange{i}) - max(e, o) + e;
    elseif i == Nproblems
        id.H{i} = 1 + max(e, o) - e:1: 1 + max(e, o) - 1;
        id.T{i} = [];
    else
        id.H{i} = 1 + max(e, o) - e:1: 1 + max(e, o) - 1;
        id.T{i} = length(alpha_subrange{i}) - max(e, o) + 1:1:length(alpha_subrange{i}) - max(e, o) + e;
    end
end

figure(1)
ax = gca;
hold(ax,'on')
plot(1:1:Nsteps, alpha_range, 'Color', 'black', 'parent', ax)
plot(horzcat(ID.h{:}),alpha_range(horzcat(ID.h{:})),'o','Color','m','markersize',16, 'parent', ax)
plot(horzcat(ID.t{:}),alpha_range(horzcat(ID.t{:})),'*','Color','c','markersize',16, 'parent', ax)
plot(horzcat(ID.H{:}),alpha_range(horzcat(ID.H{:})),'square','Color',[0,1,0.5],'markersize',30, 'parent', ax)
plot(horzcat(ID.T{:}),alpha_range(horzcat(ID.T{:})),'o','Color',[1,0,0], 'parent', ax)


end