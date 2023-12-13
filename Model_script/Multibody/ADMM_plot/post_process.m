%% post process script
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');


% Nsteps = Nsteps + 1;

X_scale = scale.x;
U_scale = scale.u;
Z_scale = scale.z;

kk = 1;
index0 = find(alpha_vec == 0);
if lap > 1
    alpha_vec2 = alpha_vec(1:index0(2)-1);
    while kk < lap
        alpha_vec2 = [alpha_vec2, alpha_vec(1:index0(2)-1) + kk];
        kk = kk + 1;
    end
    alpha_vec2 = [alpha_vec2, alpha_vec(end) + (kk-1)];
else
    alpha_vec2 = alpha_vec;
end


if Nproblems > 1
    Nsteps_0 = Nsteps/lap;
    [X_sol, U_sol, Z_sol] = denormalize_var_long(X, Nproblems, nx, nu, nz, id);
    X_sol = X_sol.*X_scale;
    U_sol = U_sol.*U_scale;
    Z_sol = Z_sol.*Z_scale;
else
    x_init = X(1:nx);
    x_other = reshape(X(nx+1:end), nu+nz+nx*(d+1), []);
    X_sol = [x_init, x_other(nu+nz+1+nx*d:end, :)].*X_scale;
    U_sol = x_other(1:nu, :).*U_scale;
    Z_sol = x_other(nu+1:nu+nz, :).*Z_scale;
    index_consensus = [];
end

time_opt = cumsum(Z_sol(3,:)); % non sono sicuro di questo calcolo



% figure
FigPlot = figure('color', 'w');
tabGrp = uitabgroup('parent', FigPlot);
labels = {'u', 'v', 'vz', 'wx', 'wy', 'r', '$\dot{l}_{11}$', '$\dot{l}_{12}$', '$\dot{l}_{21}$', '$\dot{l}_{22}$',...
    '$w_11$', '$w_12$', '$w_21$', '$w_22$', '$dgbx$', '$dgby$', '$dgbz$', '$angle_Z$', '$angle_X$', '$angle_Y$', '$l_{11}$', '$l_{12}$', '$l_{21}$', '$l_{22}$',...
    'Ta', 'Tb', 'delta',...
    'ep', 'ez', 'h', '$Fx_{11}$','$Fx_{12}$','$Fx_{21}$','$Fx_{22}$','$Fy_{11}$','$Fy_{12}$','$Fy_{21}$','$Fy_{22}$','$Fz_{11}$','$Fz_{12}$','$Fz_{21}$','$Fz_{22}$', 'trajectory'};


% plots
for i = 1:length(labels)-1
    tab = uitab(tabGrp, 'title', labels{i});
    ax = axes('parent', tab);
    hold(ax, 'on');
    if i <= nx % plot states
        plot(ax, alpha_vec2, X_sol(i, :), 'color', 'k', 'linewidth', 1.4);
        for j = 1:Nproblems
            if j == 1 
                plot(ax, alpha_vec2(ID.t{j}),  X_sol(i, ID.t{j}), 'marker', 'o', 'color', 'r', 'markersize', 15);
            elseif j == Nproblems
                plot(ax, alpha_vec2(ID.h{j}),  X_sol(i, ID.h{j}), 'marker', 'o', 'color', 'r', 'markersize', 15);
            else
                plot(ax, alpha_vec2(ID.t{j}),  X_sol(i, ID.t{j}), 'marker', 'o', 'color', 'r', 'markersize', 15);
                plot(ax, alpha_vec2(ID.h{j}),  X_sol(i, ID.h{j}), 'marker', 'o', 'color', 'r', 'markersize', 15);
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
ep_sol = Z_sol(1,:);
normals = full(pista.fun_ns(alpha_vec));
traj = pos(1:3, 2:end) + ep_sol.*normals(:,2:end);
traj = [X_sol(15:17,1),traj];
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
