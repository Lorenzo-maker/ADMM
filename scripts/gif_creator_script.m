% script for gif creation
close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
h = figure;
set(h, 'Position', get(0, 'Screensize'));
plot(middle(:,1),middle(:,2),'--',external(:,1),external(:,2),internal(:,1),internal(:,2))
title('track plot', 'FontSize', 20)
axis equal
filename = 'consensus_animation.gif';
colors = {'k', 'b', 'r', 'g','k', 'b', 'r', 'g'};
load(sprintf('%s%s',savepath, 'X_Nurb_intermediate_solutions.mat'));
linePlots = gobjects(1, Nproblems);
for i = 1:Nproblems
    linePlots(i) = line(nan,nan, 'linewidth', 1.2, 'color',colors{i});%colors{i}   
end

for n = 1:ADMM_iteration
    for j = 1:Nproblems
        
        X_current{1} = Xinter{n}{j};
        [x_unscaled] = denormalize_var_long(X_current, 0, nx, nu, nz, 0);
        Fx_opt = x_unscaled.Fx_opt;
        Fy_opt = x_unscaled.Fy_opt;
        ep_opt = x_unscaled.ep_opt;
        ef_opt = x_unscaled.ef_opt;
        alfa_opt = x_unscaled.alfa_opt;
        middle = full(pista.fun_pos(alfa_opt'))';
        x_middle = middle(:,1);
        y_middle = middle(:,2);
        normal = full(pista.fun_vh(alfa_opt'))';
        external = middle(:,1:3) + middle(:,4).*normal;
        internal = middle(:,1:3) - middle(:,4).*normal;
        vehicle = middle(:,1:3) + ep_opt.*normal;
        linePlots(j).XData = vehicle(:, 1);
        linePlots(j).YData = vehicle(:, 2);
    end
    pause(2)
    % Capture the plot as an image
    title(sprintf('track plot ADMM iteration = %i', n-1), 'FontSize', 20);
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if n == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    

end