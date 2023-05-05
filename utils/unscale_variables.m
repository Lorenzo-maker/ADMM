function [W, CONS, X, U, Z, CONS_X, CONS_U, CONS_Z] = unscale_variables(W, CONS, o, nx, nu, nz, d, dalfa)
%
% [W, CONS, X, U, Z, CONS_X, CONS_U, CONS_Z] = unscale_variables(W, CONS, nx, nu, nz, d)
%
% output: W = cell with unscaled optimization variables (states + inputs +
%               algebraic)
%         CONS = cell with unscaled consensus variables
%         X = cell with unscaled states
%         U = cell with unscaled inputs
%         Z = cell with unscaled algebraic variables
%         CONS_X = cell with unscaled consensus states
%         CONS_U = cell with unscaled consensus inputs
%         CONS_Z = cell with unscaled consensus algebraic variables
%

car_parameters_ocp;
W_scale = [ep_scale; ef_scale; d_scale; theta_scale; phi_scale; alfa_dot_scale; ep_dot_scale; ef_dot_scale; d_dot_scale; theta_dot_scale; phi_dot_scale];
U_scale = [Fx_scale; Fy_scale; Fz_scale; Mx_scale; My_scale; Mz_scale; delta_scale];
Z_scale = [Fz_scale;Fz_scale;Fz_scale;Fz_scale;Fy_scale; Fy_scale; Fx_scale; Fx_scale];

X = cell(length(W), 1);
U = cell(length(W), 1);
Z = cell(length(W), 1);

CONS_X = cell(length(W), 1);
CONS_U = cell(length(W), 1);
CONS_Z = cell(length(W), 1);

for i = 1:length(W)
    init = W{i}(1:nx).*W_scale;
    W_inter = reshape(W{i}(nx+1:end), nu+nz+nx*(d+1), []).*[U_scale; Z_scale; repmat(W_scale, d+1, 1)];
    W{i} = [init; W_inter(:)];
    
    X{i} = [init, W_inter(nu+nz+nx*d+1: end, :)]; % i-th unscaled solution for states
    U{i} = W_inter(1:nu, :);                      % i-th unscaled solution for controls    
    Z{i} = W_inter(nu+1: nu+nz, :);               % i-th unscaled solution for algebraic 
    
    if i == 1 || i == length(W); col = 1; else; col= 2; end
    CONS_inter = reshape(CONS{i}, o*(nx + nu +nz) + nx, col);
    init = CONS_inter(1:nx, :).*W_scale;
    CONS_inter = CONS_inter(nx+1:end, :).*repmat([U_scale; Z_scale; W_scale],o,1);
    CONS{i} = [init; CONS_inter];
    CONS{i} = CONS{i}(:);
    
    %CONS_X{i} = [init, CONS_inter(nu+nz+nx+1: end, :)];
    if i == 1 || i == length(W)
        temp_CONS = reshape(CONS_inter, nu+nz+nx, []);
        CONS_X{i} = [init,temp_CONS(nu+nz+1: end, :)];
    else
        temp_CONS = reshape(CONS_inter, nu+nz+nx, []);        
        %CONS_X{i} = [init(:,1),temp_CONS(nu+nz+1: end, 1:(o+1)), temp_CONS(nu+nz+1: end, (o+1)+1:end)];
        CONS_X{i} = [init(:,1),temp_CONS(nu+nz+1: end, 1:o),init(:,2), temp_CONS(nu+nz+1: end, o+1:end)];
    end
    CONS_U{i} = temp_CONS(1:nu, :);
    CONS_Z{i} = temp_CONS(nu+1: nu+nz, :);
end
