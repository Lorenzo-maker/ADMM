% check se esiste la viarabile car per calcolare twist e wrench
car_check = exist('car', 'var');
if ~car_check
    build_car;
end

for i = 1:Nproblems
    RHO_head_plots(i).XData = [RHO_head_plots(i).XData, ADMM_iteration];
    RHO_head_plots(i).YData = [RHO_head_plots(i).YData, RHO_head(i)];
    
    RHO_tail_plots(i).XData = [RHO_tail_plots(i).XData, ADMM_iteration];
    RHO_tail_plots(i).YData = [RHO_tail_plots(i).YData, RHO_tail(i)];
    
    r_head_plots(i).XData = [r_head_plots(i).XData, ADMM_iteration];
    r_head_plots(i).YData = [r_head_plots(i).YData, r_head(i)];
    
    r_tail_plots(i).XData = [r_tail_plots(i).XData, ADMM_iteration];
    r_tail_plots(i).YData = [r_tail_plots(i).YData, r_tail(i)];

    s_head_plots(i).XData = [s_head_plots(i).XData, ADMM_iteration];
    s_head_plots(i).YData = [s_head_plots(i).YData, s_head(i)];
    
    s_tail_plots(i).XData = [s_tail_plots(i).XData, ADMM_iteration];
    s_tail_plots(i).YData = [s_tail_plots(i).YData, s_tail(i)];
    
    scaled_err(i) = consErr_norm{i}./sqrt(length(Z{i}));
    scaled_change(i) = consensusChange_norm{i}./sqrt(length(Z{i}));
    
    [~, ~, states, inputs, algebr, CONS_states, CONS_inputs, CONS_algebr] = unscale_variables(X_origin, Z_1, o, nx, nu, nz, d, scale);
    
    if i == 1
        index_consensus_head = [];
        index_consensus_tail = id.t{i};
        alpha_head = [];
        alpha_tail = alpha_subrange{i}(1, index_consensus_tail);
    elseif i == Nproblems
        index_consensus_head = id.h{i};
        index_consensus_tail = [];
        alpha_head = alpha_subrange{i}(1, index_consensus_head);
        alpha_tail = [];
    else
        index_consensus_head = id.h{i};
        index_consensus_tail = id.t{i};
        alpha_head = alpha_subrange{i}(1, index_consensus_head);
        alpha_tail = alpha_subrange{i}(1, index_consensus_tail);
    end
    
    % aggiornamento delle traiettorie
    midline = full(pista.fun_pos(alpha_subrange{i}));
    normals = full(pista.fun_ns(alpha_subrange{i}));
    ep_current = states{i}(1, :);
    traj = midline(1:3, :) + ep_current.*normals;
    trajectory_plots(i).XData = traj(1, :);
    trajectory_plots(i).YData = traj(2, :);
    trajectory_plots(i).ZData = traj(3, :);
    
    % calcolo dei twist e twist al consenso
    g_gs_num = full(pista.fun_g_gs(alpha_subrange{i}));
    g_gs_num = reshape(g_gs_num, 4, 4, []);
    Jac_num = full(car.T_jac(alpha_subrange{i}));
    Jac_der_num = full(car.T_jac_der(alpha_subrange{i}));
    twist = nan(6, length(states{i}));
    for kk = 1:length(states{i})
        twist(:, kk) = full(car.axle_twist(g_gs_num(:,:,kk), Jac_num(:, kk), Jac_der_num(:, kk), states{i}(:, kk)));
    end
    
    
    g_gs_num = full(pista.fun_g_gs([alpha_head, alpha_tail]));
    g_gs_num = reshape(g_gs_num, 4, 4, []);
    Jac_num = full(car.T_jac([alpha_head, alpha_tail]));
    Jac_der_num = full(car.T_jac_der([alpha_head, alpha_tail]));
    twist_con = nan(6, 2);
    for kk = 1:length([alpha_head,alpha_tail])
        twist_con(:, kk) = full(car.axle_twist(g_gs_num(:, :, kk), Jac_num(:, kk), Jac_der_num(:, kk), CONS_states{i}(:, kk)));
    end
    
    % ciclo for per plottare gli stati
    % i = sottoproblema
    % j = stato che voglio plottare
    % gli assi partono da j+1 perchè per j=1 è associato il plot della
    % traiettoria
    % matrice states (stessa struttura anche per CONS_states):
    % [alpha_0, alpha_1, ...]
    % [   ep_0,    ep_1, ...]
    % [   ef_0,    ef_1, ...]
    % [    ...,     ..., ...]
    for j = 1:nx

        % pulisco assi
        cla(AX{i}(j+1))

        % plot di ogni stato per ogni sotto problema
        line(AX{i}(j+1), alpha_subrange{i}, states{i}(j, :));
        if i > 1 && i < Nproblems         % (i = 1 non ho head)
            line(AX{i}(j+1), alpha_head, CONS_states{i}(j, 1:o+1), 'Linestyle', 'none', 'marker', 'o', 'color', 'k');
            line(AX{i}(j+1), alpha_tail, CONS_states{i}(j, (o+1)+1:end), 'Linestyle', 'none', 'marker', 'o', 'color', 'k');
        elseif i == 1
            line(AX{i}(j+1), alpha_tail, CONS_states{i}(j, :), 'Linestyle', 'none', 'marker', 'o', 'color', 'k');
        elseif i == Nproblems
            line(AX{i}(j+1), alpha_head, CONS_states{i}(j, :), 'Linestyle', 'none', 'marker', 'o', 'color', 'k');
        end
    end
    
    % twist plot
    twist_index = [1,2,6];
    for j = nx+1:nx+3
        cla(AX{i}(j+1))
        line(AX{i}(j+1), alpha_subrange{i}, twist(twist_index(j - nx), :));
        if i > 1 && i < Nproblems %(i = 1 non ho head)
            line(AX{i}(j+1), alpha_head, twist_con(twist_index(j - nx), 1:o+1), 'Linestyle', 'none', 'marker', 'o', 'color', 'k');
            line(AX{i}(j+1), alpha_tail, twist_con(twist_index(j - nx), (o+1)+1:end), 'Linestyle', 'none', 'marker', 'o', 'color', 'k');
        elseif i == 1
            line(AX{i}(j+1), alpha_tail, twist_con(twist_index(j - nx), :), 'Linestyle', 'none', 'marker', 'o', 'color', 'k');
        elseif i == Nproblems
            line(AX{i}(j+1), alpha_head, twist_con(twist_index(j - nx), :), 'Linestyle', 'none', 'marker', 'o', 'color', 'k');
        end
    end
    
    % inputs plot
    inputs_index = [1,2,6,7];
    for j = nx+4:nx+3+4
        cla(AX{i}(j+1))
        stairs(AX{i}((j+1)), alpha_subrange{i}(1:end-1), inputs{i}(inputs_index(j-nx-3), :));
        if i > 1 && i < Nproblems %(i = 1 non ho head)
            line(AX{i}(j+1), alpha_head(1:end-1), CONS_inputs{i}(inputs_index(j-nx-3), 1:length(alpha_head)-1), 'Linestyle', 'none', 'marker', 'o', 'color', 'k');
            line(AX{i}(j+1), alpha_tail(1:end-1), CONS_inputs{i}(inputs_index(j-nx-3), length(alpha_head):end), 'Linestyle', 'none', 'marker', 'o', 'color', 'k');
        elseif i == 1
            line(AX{i}(j+1), alpha_tail(1:end-1), CONS_inputs{i}(inputs_index(j-nx-3), :), 'Linestyle', 'none', 'marker', 'o', 'color', 'k');
        elseif i == Nproblems
            line(AX{i}(j+1), alpha_head(1:end-1), CONS_inputs{i}(inputs_index(j-nx-3), :), 'Linestyle', 'none', 'marker', 'o', 'color', 'k');
        end
    end

end

% plot dei parametri di convergenza
rmax_plot.XData = [rmax_plot.XData, ADMM_iteration];
rmax_plot.YData = [rmax_plot.YData, max(scaled_err)];

smax_plot.XData = [smax_plot.XData, ADMM_iteration];
smax_plot.YData = [smax_plot.YData, max(scaled_change)];



pause(0.5)