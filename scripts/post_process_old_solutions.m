%% load solutions and post process
clc; clear all;
reults_name = '1500_1_0';

addpath(genpath('Casadi'));
addpath(genpath('ADMM_functions'));
addpath(genpath('Classes'));
addpath(genpath('Track'));
addpath(genpath('utils'));
addpath(genpath('Data'));

%% loading results
path = sprintf('Results\\%s', reults_name);

load(sprintf('%s\\%s', path, 'opti_var_scaled.mat'));
load(sprintf('%s\\%s', path, 'X_sol.mat'));
load(sprintf('%s\\%s', path, 'Z_sol.mat'));
load(sprintf('%s\\%s', path, 'U_sol.mat'));
load(sprintf('%s\\%s', path, 'problem_structure.mat'));

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');

%% extracting necessary variables



ADMM_batch_settings;

alpha_subrange = problem_structure.alpha_subrange;
overlap = problem_structure.overlap;
alfa_end = problem_structure.alpha_end;
alpha_vec = problem_structure.alpha_vec;
dalfa = alpha_vec(2)-alpha_vec(1);
Nsteps = problem_structure.Nsteps;
Nproblems = problem_structure.Nproblems;
car_parameters_ocp;
load('Data\pista.mat');
build_car;

%% plots
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
alpha_vec_colloc = linspace(0, alfa_end, ((d+1))*(Nsteps-1)+1);
if Nproblems>1
    [index_consensus, index_end] = alphaSubToIndex(alpha_subrange, overlap, Nproblems, alfa_end, Nsteps+1);
    
    [X_sol, U_sol, Z_sol] = denormalize_var_long(X, Nproblems, nx, nu, nz,  overlap, dalfa);
    alpha_vec = linspace(0, alfa_end, Nsteps+1);
    X_sol = [alpha_vec; X_sol].*X_scale;
    U_sol = U_sol.*U_scale;
    Z_sol = Z_sol.*Z_scale;
else
    alpha_vec = linspace(0, alfa_end, Nsteps);
    x_init = X(1:nx);
    x_other = reshape(X(nx+1:end), nu+nz+nx*(d+1), []);
    X_sol = [alpha_vec; x_init, x_other(nu+nz+1+nx*d:end, :)].*X_scale;
    other_states = reshape(x_other(nu+nz+1:end, :), nx, (Nsteps-1)*(d+1) ).*X_scale(2:end);
    X_sol_colloc = [x_init.*X_scale(2:end), other_states];
    X_sol_colloc = [alpha_vec_colloc.*X_scale(1); X_sol_colloc];
    U_sol = x_other(1:nu, :).*U_scale;
    Z_sol = x_other(nu+1:nu+nz, :).*Z_scale;
    index_consensus = [];
end


time_opt = sum((1./X_sol(7,2:end)).*dalfa); % non sono sicuro di questo calcolo

% compute twist
g_gs_num = reshape(full(pista.fun_g_gs(alpha_vec)), 4, 4, []);
g_gs_num_colloc = reshape(full(pista.fun_g_gs(alpha_vec_colloc)), 4, 4, []);
J_num = full(car.T_jac(alpha_vec));
J_dot_num = full(car.T_jac_der(alpha_vec));
twist = nan(6, size(g_gs_num, 3));
for i = 1:size(g_gs_num, 3)
    twist(:, i) = full(car.axle_twist(g_gs_num(:,:,i), J_num(:, i), J_dot_num(:, i), X_sol(2:end, i)));
end

% figure
FigPLot = figure('color', 'w');
tabGrp = uitabgroup('parent', FigPLot);
labels = {'ep', 'ef', 'd', 'theta', 'phi', 'alpha dot', 'ep dot', 'ef dot', 'd dot', 'theta dot', 'phi dot',...
    'u', 'v', 'w', 'omega_x', 'omega_y', 'r', 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz', 'delta', 'trajectory'};


% plots
for i = 1:length(labels)-1
    tab = uitab(tabGrp, 'title', labels{i});
    ax = axes('parent', tab);
    hold(ax, 'on');
    if i <= nx % plot states
        plot(ax, alpha_vec, X_sol(i+1, :), 'color', 'k', 'linewidth', 1.4);
        for j = 1:length(index_consensus)
            plot(ax, alpha_vec(index_consensus(j)),  X_sol(i+1, index_consensus(j)), 'marker', 'o', 'color', 'r', 'markersize', 15);
        end
    elseif i > nx && i<= nx+6 % plot twist
        ii = i-nx;
        plot(ax, alpha_vec, twist(ii, :), 'color', 'k', 'linewidth', 1.4);
        for j = 1:length(index_consensus)
            plot(ax, alpha_vec(index_consensus(j)),  twist(ii, index_consensus(j)), 'marker', 'o', 'color', 'r', 'markersize', 15);
        end
    elseif i > nx+6 && i <= nx+6+nu
        ii = i-nx-6;
        plot(ax, alpha_vec(2:end), U_sol(ii, :), 'color', 'k', 'linewidth', 1.4);
        for j = 1:length(index_consensus)
            plot(ax, alpha_vec(index_consensus(j)),  U_sol(ii, index_consensus(j)-1), 'marker', 'o', 'color', 'r', 'markersize', 15);
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
ep_sol_colloc = X_sol_colloc(2,:);
n = permute(g_gs_num(1:3, 2, :), [1 3 2]);
n_colloc = permute(g_gs_num_colloc(1:3, 2, :), [1 3 2]);
w = permute(g_gs_num(1:3, 3, :), [1 3 2]);
w_colloc = permute(g_gs_num_colloc(1:3, 3, :), [1 3 2]);
pos_colloc = full(pista.fun_pos(alpha_vec_colloc));
traj_colloc = pos_colloc(1:3, :) + ep_sol_colloc.*n_colloc + 0.001.*w_colloc;
traj = pos(1:3, :) + ep_sol.*n + 0.001.*w;
line(traj(1, :), traj(2, :), traj(3, :), 'color', 'k', 'linewidth', 2, 'parent', ax)
set(ax, 'Fontsize', 18)
xlabel('x (m)')
ylabel('y (m)')
set(ax, 'DataAspectRatio', [1 1 1])
view(0,90)
