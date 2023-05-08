%% post process script
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');


% Nsteps = Nsteps + 1;
alpha_vec2 = [];
kk = 0;
if Nproblems > 1
    %[index_consensus, index_end] = alphaSubToIndex(alpha_subrange, overlap, Nproblems, alfa_end, Nsteps+1);
    Nsteps_0 = Nsteps/lap;
    [X_sol, U_sol, Z_sol] = denormalize_var_long(X, Nproblems, nx, nu, nz, id);
    if lap > 1
        alpha_vec_lap = linspace(0, alfa_end, (Nsteps_0 + 1));
        alpha_vec = repmat(alpha_vec_lap(1:end-1), 1, lap-1);
        alpha_vec = [alpha_vec, alpha_vec_lap];
        while kk < lap
            alpha_vec2 = [alpha_vec2, alpha_vec_lap(1:end-1) + kk];
            kk = kk + 1;
        end
        alpha_vec2(end+1) = lap;
    else
        alpha_vec = linspace(0, alfa_end, Nsteps+1);
        alpha_vec2 = alpha_vec;
    end
    %alpha_vec = linspace(0, alfa_end, Nsteps+1);
    X_sol = [alpha_vec; X_sol].*[1;X_scale];
    U_sol = U_sol.*U_scale;
    Z_sol = Z_sol.*Z_scale;
else
    if lap > 1
        alpha_vec_lap = linspace(0, alfa_end, (Nsteps_0 + 1));
        alpha_vec = repmat(alpha_vec_lap(1:end-1), 1, lap-1);
        alpha_vec = [alpha_vec, alpha_vec_lap];
        while kk < lap
            alpha_vec2 = [alpha_vec2, alpha_vec_lap(1:end-1) + kk];
            kk = kk + 1;
        end
        alpha_vec2(end+1) = lap;
    else
        alpha_vec = linspace(0, alfa_end, Nsteps);
        alpha_vec2 = alpha_vec;
    end
    %alpha_vec = linspace(0, alfa_end, Nsteps+1);
    alpha_vec = linspace(0, alfa_end, Nsteps);
    x_init = X(1:nx);
    x_other = reshape(X(nx+1:end), nu+nz+nx*(d+1), []);
    if alpha_numeric
        X_sol = [alpha_vec; x_init, x_other(nu+nz+1+nx*d:end, :)].*[1;X_scale];
        %other_states = reshape(x_other(nu+nz+1:end, :), nx, (Nsteps-1)*(d+1) ).*X_scale(2:end);
    else
        X_sol = [x_init, x_other(nu+nz+1+nx*d:end, :)].*X_scale;
    end
    U_sol = x_other(1:nu, :).*U_scale;
    Z_sol = x_other(nu+1:nu+nz, :).*Z_scale;
    index_consensus = [];
end


time_opt = cumsum((1./X_sol(7,2:end)).*dalfa); % non sono sicuro di questo calcolo



% figure
FigPlot = figure('color', 'w');
tabGrp = uitabgroup('parent', FigPlot);
labels = {'u', 'v', 'r', 'x', 'y', 'psi', 'Fx', 'rp', 'Fz', 'Fy', 'P', 'ep', 'h', 'trajectory'};


% plots
for i = 1:length(labels)-1
    tab = uitab(tabGrp, 'title', labels{i});
    ax = axes('parent', tab);
    hold(ax, 'on');
    if i <= nx % plot states
        plot(ax, alpha_vec2, X_sol(i+1, :), 'color', 'k', 'linewidth', 1.4);
        for j = 1:Nproblems
            if j == 1 
                plot(ax, alpha_vec2(ID.t{j}),  X_sol(i+1, ID.t{j}), 'marker', 'o', 'color', 'r', 'markersize', 15);
            elseif j == Nproblems
                plot(ax, alpha_vec2(ID.h{j}),  X_sol(i+1, ID.h{j}), 'marker', 'o', 'color', 'r', 'markersize', 15);
            else
                plot(ax, alpha_vec2(ID.t{j}),  X_sol(i+1, ID.t{j}), 'marker', 'o', 'color', 'r', 'markersize', 15);
                plot(ax, alpha_vec2(ID.h{j}),  X_sol(i+1, ID.h{j}), 'marker', 'o', 'color', 'r', 'markersize', 15);
            end
        end
    elseif i > nx && i <= nx+nu
        ii = i-nx;
        plot(ax, alpha_vec2(1:end-1), U_sol(ii, :), 'color', 'k', 'linewidth', 1.4);
        for j = 1:Nproblems
            if j == 1 
                plot(ax, alpha_vec2(ID.t{j}(1:end-1)),  U_sol(ii, ID.t{j}(1:end-1)), 'marker', 'o', 'color', 'r', 'markersize', 15);
            elseif j == Nproblems
                plot(ax, alpha_vec2(ID.h{j}(1:end-1)),  U_sol(ii, ID.h{j}(1:end-1)), 'marker', 'o', 'color', 'r', 'markersize', 15);
            else
                plot(ax, alpha_vec2(ID.t{j}(1:end-1)),  U_sol(ii, ID.t{j}(1:end-1)), 'marker', 'o', 'color', 'r', 'markersize', 15);
                plot(ax, alpha_vec2(ID.h{j}(1:end-1)),  U_sol(ii, ID.h{j}(1:end-1)), 'marker', 'o', 'color', 'r', 'markersize', 15);
            end
        end
    elseif i > nx+nu && i <= nx+nu+nz
        ii = i-nx-nu;
        plot(ax, alpha_vec2(1:end-1), Z_sol(ii, :), 'color', 'k', 'linewidth', 1.4);
        for j = 1:Nproblems
            if j == 1 
                plot(ax, alpha_vec2(ID.t{j}(1:end-1)),  Z_sol(ii, ID.t{j}(1:end-1)), 'marker', 'o', 'color', 'r', 'markersize', 15);
            elseif j == Nproblems
                plot(ax, alpha_vec2(ID.h{j}(1:end-1)),  Z_sol(ii, ID.h{j}(1:end-1)), 'marker', 'o', 'color', 'r', 'markersize', 15);
            else
                plot(ax, alpha_vec2(ID.t{j}(1:end-1)),  Z_sol(ii, ID.t{j}(1:end-1)), 'marker', 'o', 'color', 'r', 'markersize', 15);
                plot(ax, alpha_vec2(ID.h{j}(1:end-1)),  Z_sol(ii, ID.h{j}(1:end-1)), 'marker', 'o', 'color', 'r', 'markersize', 15);
            end
        end
    end    
    set(ax, 'Fontsize', 18)
    xlabel('$\alpha$')
    ylabel(labels{i})
end

% trajectory plot
tab = uitab(tabGrp, 'title', 'Trajectory');
ax = axes('parent', tab);
hold (ax, 'on')
pos = full(pista.fun_pos(alpha_vec));
pista.full_plot(alpha_vec, 'dotted', false, 'parent', ax)
ep_sol = [Z_sol(4,:),Z_sol(4,end)];
n = permute(g_gs_num(1:3, 2, :), [1 3 2]);
w = permute(g_gs_num(1:3, 3, :), [1 3 2]);
traj = pos(1:3, :) + ep_sol.*n + 0.001.*w;
line(traj(1, :), traj(2, :), traj(3, :), 'color', 'k', 'linewidth', 2, 'parent', ax)
set(ax, 'Fontsize', 18)
xlabel('x (m)')
ylabel('y (m)')
set(ax, 'DataAspectRatio', [1 1 1])
view(0,90)

%
% figure(2)
% plot(alpha_vec(1:end-1), Z_sol(1,:), alpha_vec(1:end-1), Z_sol(2,:))
% 
% figure(3)
% plot(alpha_vec(1:end-1), twist(1,1:end-1).*U_sol(1,:))
% hold on
% line([0,1],47.1e3*[1,1],'Color','red','LineStyle','--', 'Linewidth', 1.2)
