%% post process script
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');


% Nsteps = Nsteps + 1;
alpha_vec2 = [];
kk = 0;
if Nproblems > 1
    [index_consensus, index_end] = alphaSubToIndex(alpha_subrange, overlap, Nproblems, alfa_end, Nsteps+1);
    Nsteps_0 = Nsteps/lap;
    [X_sol, U_sol, Z_sol] = denormalize_var_long(X, Nproblems, nx, nu, nz,  overlap, dalfa);
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
    X_sol = [alpha_vec; X_sol].*X_scale;
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
    X_sol = [alpha_vec; x_init, x_other(nu+nz+1+nx*d:end, :)].*X_scale;
    other_states = reshape(x_other(nu+nz+1:end, :), nx, (Nsteps-1)*(d+1) ).*X_scale(2:end);
    U_sol = x_other(1:nu, :).*U_scale;
    Z_sol = x_other(nu+1:nu+nz, :).*Z_scale;
    index_consensus = [];
end


time_opt = cumsum((1./X_sol(7,2:end)).*dalfa); % non sono sicuro di questo calcolo

% compute twist
g_gs_num = reshape(full(pista.fun_g_gs(alpha_vec)), 4, 4, []);
J_num = full(car.T_jac(alpha_vec));
J_dot_num = full(car.T_jac_der(alpha_vec));
twist = nan(6, size(g_gs_num, 3));
for i = 1:size(g_gs_num, 3)
    twist(:, i) = full(car.axle_twist(g_gs_num(:,:,i), J_num(:, i), J_dot_num(:, i), X_sol(2:end, i)));
end

% figure
FigPlot = figure('color', 'w');
tabGrp = uitabgroup('parent', FigPlot);
labels = {'ep', 'ef', 'd', 'theta', 'phi', 'alpha dot', 'ep dot', 'ef dot', 'd dot', 'theta dot', 'phi dot',...
    'u', 'v', 'w', 'omega_x', 'omega_y', 'r', 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz', 'delta', 'trajectory'};


% plots
for i = 1:length(labels)-1
    tab = uitab(tabGrp, 'title', labels{i});
    ax = axes('parent', tab);
    hold(ax, 'on');
    if i <= nx % plot states
        plot(ax, alpha_vec2, X_sol(i+1, :), 'color', 'k', 'linewidth', 1.4);
        for j = 1:length(index_consensus)
            plot(ax, alpha_vec2(index_consensus(j)),  X_sol(i+1, index_consensus(j)), 'marker', 'o', 'color', 'r', 'markersize', 15);
        end
    elseif i > nx && i<= nx+6 % plot twist
        ii = i-nx;
        plot(ax, alpha_vec2, twist(ii, :), 'color', 'k', 'linewidth', 1.4);
        for j = 1:length(index_consensus)
            plot(ax, alpha_vec2(index_consensus(j)),  twist(ii, index_consensus(j)), 'marker', 'o', 'color', 'r', 'markersize', 15);
        end
    elseif i > nx+6 && i <= nx+6+nu
        ii = i-nx-6;
        plot(ax, alpha_vec2(2:end), U_sol(ii, :), 'color', 'k', 'linewidth', 1.4);
        for j = 1:length(index_consensus)
            plot(ax, alpha_vec2(index_consensus(j)),  U_sol(ii, index_consensus(j)-1), 'marker', 'o', 'color', 'r', 'markersize', 15);
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
ep_sol = X_sol(2,:);
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