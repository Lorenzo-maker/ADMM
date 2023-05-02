function out = update_var_acc(v_bar, v_bar_prev, v_acc_all, nprob, iter, options)

arguments
   v_bar 
   v_bar_prev
   v_acc_all
   nprob
   iter
   options.method = 'Stationary'; 
   options.alfa = 0.25;
end
beta(1) = 1;
for i = 1:iter+1
    beta = [beta;(1 + sqrt(1+ 4*beta(i)^2))/2];
end
    

switch options.method
    case 'ADMM'
        out = v_bar;
    case 'Stationary'
        alfa = options.alfa;
        out = v_bar + alfa*(v_bar - v_bar_prev);
    case 'Nesterov' % norm 1
        
        alfa = (beta(iter + 1) - 1)/beta(iter + 2);
        out = v_bar + alfa*(v_bar - v_bar_prev);   
    case 'Automatic' % norm 1   
        if iter < 2
            alfa = 0;
            out = v_bar + alfa*(v_bar - v_bar_prev);
        else
            gk_1 = (v_bar - vertcat(v_acc_all{iter}{nprob}));
            gk_2 = (v_bar_prev - vertcat(v_acc_all{iter-1}{nprob}));
            alfa = gk_1'*gk_2/(gk_2'*gk_2);
            if alfa < 0
                alfa = 0;
            elseif alfa > 1
                    alfa = 1;
            end
            out = v_bar + alfa*(v_bar - v_bar_prev);
        end
    otherwise
        error('select appropriate update method')
end

end