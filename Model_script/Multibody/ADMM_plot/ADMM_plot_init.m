%% plots init

% latex interpreter appearence
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex'); 
set(0,'defaultTextInterpreter','latex');
% figure
F = figure('color', 'w');
tabG = uitabgroup(F);

tabs = gobjects(Nproblems+1, 1);
subtabtitle = {'convergence', 'trajectory', 'u', 'v', 'vz', 'wx', 'wy', 'r', '$\dot{l}_{11}$', '$\dot{l}_{12}$', '$\dot{l}_{21}$', '$\dot{l}_{22}$',...
    '$w_11$', '$w_12$', '$w_21$', '$w_22$', '$dgbx$', '$dgby$', '$dgbz$', '$angle_Z$', '$angle_X$', '$angle_Y$', '$l_{11}$', '$l_{12}$', '$l_{21}$', '$l_{22}$',...
    'Ta', 'Tb', 'delta',...
    'ep', 'ez', 'h', '$Fx_{11}$','$Fx_{12}$','$Fx_{21}$','$Fx_{22}$','$Fy_{11}$','$Fy_{12}$','$Fy_{21}$','$Fy_{22}$','$Fz_{11}$','$Fz_{12}$','$Fz_{21}$','$Fz_{22}$'};

for i = 1:Nproblems+1
    
    if i == 1
        tabtitle = sprintf('main'); 
    else
        tabtitle = sprintf('sub-problem %i', i-1);
    end
    tabs(i) = uitab(tabG, 'title', tabtitle);
    
    if i > 1
        tabSubGrp = uitabgroup(tabs(i));
        for j = 1:length(subtabtitle)
            subtabs{i-1}(j) = uitab(tabSubGrp, 'title', subtabtitle{j});
        end
    end
end

%% axis creation for ADMM convergence plots
main_ax1 = subplot( 2, 1, 1); hold on; box on; set(main_ax1, 'yscale', 'log'); set(main_ax1, 'parent', tabs(1)); % max primal res
main_ax2 = subplot( 2, 1, 2); hold on; box on; set(main_ax2, 'yscale', 'log'); set(main_ax2, 'parent', tabs(1)); % max dual res

title(main_ax1, 'max consensus error');
xlabel(main_ax1, 'ADMM iteration');
ylabel(main_ax1, '$\frac{||\bf{r}_\textrm{max}||_2}{\sqrt{n}}$');

title(main_ax2, 'max consensus change');
xlabel(main_ax2, 'ADMM iteration');
ylabel(main_ax2, '$\frac{||\bf{s}_\textrm{max}||_2}{\sqrt{n}}$');

set(get(tabs(1), 'children'), 'Fontsize', 16)

rmax_plot = plot(main_ax1, 0, 0, 'linewidth', 1.5);
rmax_plot.XData = nan;
rmax_plot.YData = nan;
smax_plot = plot(main_ax2, 0, 0, 'linewidth', 1.5);
smax_plot.XData = nan;
smax_plot.YData = nan;
    
RHO_tail_plots = gobjects(Nproblems, 1);
RHO_head_plots = gobjects(Nproblems, 1);
r_head_plots = gobjects(Nproblems, 1);
r_tail_plots = gobjects(Nproblems, 1);
s_head_plots = gobjects(Nproblems, 1);
s_tail_plots = gobjects(Nproblems, 1);

problem_axis = gobjects(Nproblems, 6);
for i = 1:Nproblems
    
    for k = 1:6
        % 1 penalty head axis
        % 2 penalty tail axis
        % 3 primal res head axis
        % 4 primal res tail axis
        % 5 dual res head axis
        % 6 dual res tail axis
        problem_axis(i, k) = subplot( 3, 2, k, 'parent', subtabs{i}(1)); hold(problem_axis(i, k), 'on'); box on; set(problem_axis(i, k), 'yscale', 'log');
        xlabel(problem_axis(i, k), 'ADMM iteration');
    end
    
    title(problem_axis(i, 1), 'head penalty');
    ylabel(problem_axis(i, 1), '$\rho_h$');

    title(problem_axis(i, 2), 'tail penalty');
    ylabel(problem_axis(i, 2), '$\rho_t$');

    title(problem_axis(i, 3), 'primal residual head');
    ylabel(problem_axis(i, 3), '$||\bf{r}_\textrm{h}||$');

    title(problem_axis(i, 4), 'primal residual tail');
    ylabel(problem_axis(i, 4), '$||\bf{r}_\textrm{t}||$');

    title(problem_axis(i, 5), 'dual residual head');
    ylabel(problem_axis(i, 5), '$||\bf{s}_\textrm{h}||$');

    title(problem_axis(i, 6), 'dual residual tail');
    ylabel(problem_axis(i, 6), '$||\bf{s}_\textrm{t}||$');

    
    RHO_head_plots(i) = plot(problem_axis(i, 1), 0, 0, 'linewidth', 1.5);
    RHO_head_plots(i).XData = nan;
    RHO_head_plots(i).YData = nan;
    
    RHO_tail_plots(i) = plot(problem_axis(i, 2), 0, 0, 'linewidth', 1.5);
    RHO_tail_plots(i).XData = nan;
    RHO_tail_plots(i).YData = nan;
    
    r_head_plots(i) = plot(problem_axis(i, 3), 0, 0, 'linewidth', 1.5);
    r_head_plots(i).XData = nan;
    r_head_plots(i).YData = nan;
    
    r_tail_plots(i) = plot(problem_axis(i, 4), 0, 0, 'linewidth', 1.5);
    r_tail_plots(i).XData = nan;
    r_tail_plots(i).YData = nan;
    
    s_head_plots(i) = plot(problem_axis(i, 5), 0, 0, 'linewidth', 1.5);
    s_head_plots(i).XData = nan;
    s_head_plots(i).YData = nan;
    
    s_tail_plots(i) = plot(problem_axis(i, 6), 0, 0, 'linewidth', 1.5);
    s_tail_plots(i).XData = nan;
    s_tail_plots(i).YData = nan;
    
    set(get(subtabs{i}(1), 'children'), 'Fontsize', 16)
end

%% axis creation for intermediate solutions plot
AX = cell(Nproblems, 1);
trajectory_plots = gobjects(Nproblems, 1);
for i = 1:Nproblems
    AX{i} = gobjects(length(subtabtitle)-1, 1);
    for j = 2:length(subtabtitle)

        AX{i}(j-1) = axes(subtabs{i}(j), 'box', 'on');
        if j == 2
            axis equal; hold on;
            % init track sector plot
            pista.full_plot(alpha_subrange{i}, 'patch', 0, 'parent', AX{i}(j-1));
            trajectory_plots(i) = line(AX{i}(j-1), nan, nan, nan, 'linewidth', 1.1, 'color', 'k');
        else
            xlabel(AX{i}(j-1), '$\alpha$');
        end
        set(0,'defaultTextInterpreter','latex');
        title(AX{i}(j-1), subtabtitle{j});
    end

end


