function Z = update_dual_var(X_tail, X_head, RHO_tail, RHO_head, Y_tail, Y_head, options)
%
% solution of the augmented lagrangian w.r.t. the dual variable Z,
% considering different penalties coefficients
%
arguments
    X_tail
    X_head
    RHO_tail
    RHO_head
    Y_tail
    Y_head
    options.method = 'Augmented-Lagrangian';
end
switch options.method
    case 'Augmented-Lagrangian'
        Z = (Y_tail + Y_head + RHO_tail.*X_tail + RHO_head.*X_head)./(RHO_tail + RHO_head);
    case 'exact-penalty' % norm 1
        Z = (abs(RHO_tail.*X_tail) > abs(RHO_head.*X_head)).*X_head + (abs(RHO_tail.*X_tail) < abs(RHO_head.*X_head)).*X_tail;
        if  RHO_tail == RHO_head % infinite solutions
            Z = (X_tail + X_head)./2;
        end
    otherwise
        error('select appropriate update method')
end
end

