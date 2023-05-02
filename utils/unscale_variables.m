function [W, CONS, X, U, Z, CONS_X, CONS_U, CONS_Z] = unscale_variables(W, CONS, nx, nu, nz, d, dalfa)
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
    
    X{i} = [init, W_inter(nu+nz+nx*d+1: end, :)];
    U{i} = W_inter(1:nu, :);
    Z{i} = W_inter(nu+1: nu+nz, :);
    
    if i == 1 || i == length(W); col = 1; else; col= 2; end
    CONS_inter = reshape(CONS{i}, nx*2 + nu +nz, col);
    init = CONS_inter(1:nx, :).*W_scale;
    CONS_inter = CONS_inter(nx+1:end, :).*[U_scale; Z_scale; W_scale];
    CONS{i} = [init; CONS_inter];
    CONS{i} = CONS{i}(:);
    
    %CONS_X{i} = [init, CONS_inter(nu+nz+nx+1: end, :)];
    if i == 1 || i == length(W)
        CONS_X{i} = [init,CONS_inter(nu+nz+1: end, :)];
    else
        CONS_X{i} = [init(:,1),CONS_inter(nu+nz+1: end, 1),init(:,2), CONS_inter(nu+nz+1: end, 2)];
    end
    CONS_U{i} = CONS_inter(1:nu, :);
    CONS_Z{i} = CONS_inter(nu+1: nu+nz, :);
end
