%% post process script
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
clear all
load('Results/1500_4_81/ws_4_11')
[index_consensus, index_end, index_local] = alphaSubToIndex(alpha_subrange, overlap, Nproblems, alfa_end, Nsteps+1);
close all
filename = 'consensus_animation.gif';
kkp = [5,6,7,8,9,13,18,20,21];
kkn = [1,2,3,4,5,6,8,11,13];
count = 1;


%%%%%%%%%%%%%%%%%
%2*ep track
%%%%%%%%%%%%%%%%%

nn = 9
for kk = kkp
    X = Xinter{kk};
    Z = Zinter{kk};
    data = cell(Nproblems,1);
    for i = 1:Nproblems
        ZZ = reshape(Z{i}, nx+nu+nz+nx, []);
        ZZ = ZZ.*[X_scale(2:end);U_scale;Z_scale;X_scale(2:end)];
        XX = reshape(X{i}(nx+1:end),nu+nz+nx,[]).*[U_scale;Z_scale;X_scale(2:end)];
        XX0 = X{i}(1:nx).*X_scale(2:end);
        data{i}.Z = ZZ;
        data{i}.sol = XX;
        data{i}.sol_0 = XX0;
        data{i}.alpha = alpha_subrange{i};    
        if i < Nproblems
            data{i}.index_consensus = index_consensus(i);
        end
        data{i}.g_gs_num = reshape(full(pista.fun_g_gs(alpha_subrange{i})), 4, 4, []);
        data{i}.J_num = full(car.T_jac(alpha_subrange{i}));
        data{i}.J_dot_num = full(car.T_jac_der(alpha_subrange{i}));
        data{i}.twist = nan(6, size(data{i}.g_gs_num, 3));
        for j = 1:size(data{i}.g_gs_num, 3)
            X_sol = [XX0,XX(nu+nz+1:end,:)];
            data{i}.twist(:, j) = full(car.axle_twist(data{i}.g_gs_num(:,:,j), data{i}.J_num(:, j), data{i}.J_dot_num(:, j), X_sol(:, j)));
        end
        data{i}.pos = full(pista.fun_pos(alpha_subrange{i}));
        ep_sol = 2*X_sol(1,:);
        n = permute(data{i}.g_gs_num(1:3, 2, :), [1 3 2]);
        w = permute(data{i}.g_gs_num(1:3, 3, :), [1 3 2]);
        data{i}.traj = data{i}.pos(1:3, :) + ep_sol.*n + 0.001.*w;
        % Computed consensus index for each subproblems
        if i > 1 && i < Nproblems
           alpha_head = alpha_subrange{i-1}(index_local(i-1));
           alpha_tail = alpha_subrange{i}(index_local(i));
           data{i}.index_head = find(alpha_subrange{i} == alpha_head);
           data{i}.index_tail = find(alpha_subrange{i} == alpha_tail);
        elseif i == 1
           data{i}.index_tail = index_local(i);
        elseif i == Nproblems
           data{i}.index_head = find(alpha_subrange{i} == alpha_tail);
        end            
        % Compute consensus on trajectory
        if i < Nproblems
            pos = data{i}.pos(data{i}.index_tail:data{i}.index_tail+1);
            epZ = 2*[data{i}.Z(1,end), data{i}.Z(nx+nu+nz+1,end)];
            data{i}.trajZ = data{i}.pos(1:3, data{i}.index_tail:data{i}.index_tail+1) + epZ.*n(:,data{i}.index_tail:data{i}.index_tail+1) + 0.001.*w(:,data{i}.index_tail:data{i}.index_tail+1);
        end




    end

    % figure(1)
    % ax = gca;
    % hold(ax, 'on')
    % plot(data{1}.alpha(1:index_local(1)), [data{1}.sol_0(1), data{1}.sol(nu+nz+1,1:index_local(1)-1)],'marker','x','parent',ax)
    % plot(data{1}.alpha(index_consensus(1):index_consensus(1)+1), data{1}.sol(nu+nz+1,index_consensus(1)-1:index_consensus(1)),'marker', 'o')
    % plot(data{1}.alpha(index_local(1)+1:end), data{1}.sol(nu+nz+1,index_local(1):end),'marker','x','parent',ax)
    % plot(data{2}.alpha, [data{2}.sol_0(2), data{2}.sol(nu+nz+1,:)],'parent',ax)
    h = figure(2);
    set(h, 'Position', get(0, 'Screensize'));
    count = count + 1;
    ax2 = gca;
    hold (ax2, 'on')
    pista.full_plot_admm(linspace(0,1,1501), 'dotted', false, 'parent', ax2)
    view(0, 90) %view(74.157725341694245,79.113322170517506)
%     camproj(ax2, 'perspective')
%     camva(ax2, 3)
%     camtarget(ax2, [-2030 2750 -10])
    xlim([-2090 -1980])%xlim([-2100 -1980])
    ylim([2690 2810])%ylim([2690 2810])
    %zlim([-200 0])
    set(ax2,'DataAspectRatio',[1 1 1])
    ax3 = ax2;
    axis('off')
    %axis('equal')
    grid on
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(0,'defaultTextInterpreter','latex');
    for k = 1:Nproblems
        if k > 1 && k < Nproblems
            if k == 2
                visible = 'on';
            else
                visible = 'off';
            end
            plot(data{k}.traj(1,3:data{k}.index_head+1),data{k}.traj(2,3:data{k}.index_head+1),'color',[0.45,0.95,0.45],'marker','o','MarkerFaceColor',[0.45,0.95,0.45],'parent', ax3, 'Linewidth', 8,'HandleVisibility',visible,'DisplayName','traj$_{2,H}$')
            plot(data{k}.traj(1,data{k}.index_head+1:data{k}.index_tail),data{k}.traj(2,data{k}.index_head+1:data{k}.index_tail),'black','marker','o','MarkerFaceColor','black','parent', ax3, 'Linewidth', 2,'HandleVisibility',visible,'DisplayName','traj$_{2}$')
            plot(data{k}.traj(1,data{k}.index_head:data{k}.index_head+1),data{k}.traj(2,data{k}.index_head:data{k}.index_head+1),'m','marker','o','MarkerFaceColor','m','parent', ax3, 'Linewidth', 4,'HandleVisibility',visible,'DisplayName','traj$_{2,h}$')            
            plot(data{k}.traj(1,data{k}.index_tail:end),data{k}.traj(2,data{k}.index_tail:end),'r','marker','o','MarkerFaceColor','r','Linewidth', 2,'HandleVisibility','off')
            plot(data{k}.traj(1,data{k}.index_tail:data{k}.index_tail+1),data{k}.traj(2,data{k}.index_tail:data{k}.index_tail+1),'c','marker','o','MarkerFaceColor','c', 'Linewidth', 2,'parent', ax3,'HandleVisibility','off')                       
            plot(data{k}.trajZ(1,:),data{k}.trajZ(2,:),'color',[0.9290 0.6940 0.1250],'marker','diamond','markersize',12,'MarkerFaceColor',[0.9290 0.6940 0.1250],'parent', ax3,'HandleVisibility','off')
        elseif k == 1
            plot(data{k}.trajZ(1,:),data{k}.trajZ(2,:),'color',[0.9290 0.6940 0.1250],'marker','diamond','markersize',12,'MarkerFaceColor',[0.9290 0.6940 0.1250],'parent', ax3, 'Linewidth', 2,'DisplayName','traj$_{z_{1}}$')
            plot(data{k}.traj(1,1:data{k}.index_tail),data{k}.traj(2,1:data{k}.index_tail),'black','marker','square','MarkerFaceColor','black','parent', ax3, 'Linewidth', 2,'DisplayName','traj$_{1}$')
            plot(data{k}.traj(1,data{k}.index_tail:end-2),data{k}.traj(2,data{k}.index_tail:end-2),'color',[1,0.3,0.45],'marker','o','MarkerFaceColor',[0.45,0.95,0.45],'parent', ax3, 'Linewidth', 8,'DisplayName','traj$_{1,T}$')
            plot(data{k}.traj(1,data{k}.index_tail:data{k}.index_tail+1),data{k}.traj(2,data{k}.index_tail:data{k}.index_tail+1),'c','marker','o','MarkerFaceColor','c', 'Linewidth', 4,'parent', ax3,'DisplayName','traj$_{1,t}$')
        elseif k == Nproblems
            plot(data{k}.traj(1,1:data{k}.index_head),data{k}.traj(2,1:data{k}.index_head),'g','marker','o','MarkerFaceColor','g','parent', ax3,'HandleVisibility','off')
            plot(data{k}.traj(1,data{k}.index_head:end),data{k}.traj(2,data{k}.index_head:end),'black','marker','o','MarkerFaceColor','black','parent', ax3,'HandleVisibility','off') 
            plot(data{k}.traj(1,data{k}.index_head:data{k}.index_head+1),data{k}.traj(2,data{k}.index_head:data{k}.index_head+1),'m','marker','o','MarkerFaceColor','m','parent', ax3,'HandleVisibility','off')
        end
    end
    legend('Location','none', 'Position', [0.6,0.37,0.19,0.59],'FontSize', 34)
    text(-2080, 2722, sprintf('$k = %i$', kkn(count - 1)), 'FontSize', 34)
    pause(2)       
%     title(ax3,sprintf('track plot ADMM iteration = %i', kkn(count-1)), 'FontSize', 20);
%     text(-2080, 2700,-30,'\leftarrow sin(\pi)')
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if kk == 5
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    if kk < ADMM_iteration
     cla(ax3)
    end
    kk
 end