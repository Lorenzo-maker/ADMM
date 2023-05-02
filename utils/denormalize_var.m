function [x_unscaled] = denormalize_var(X, Nproblems, nx, nu, nz, elemOverlap)
%
% x_unscaled = denormalize_var(x_scaled)
%

car_parameters;

% concatenate solutions
sol = [];
sol = cat(1, sol, X{1});
for i = 2:Nproblems
    
    sol = cat(1, sol, X{i}(1 + elemOverlap: end));
    
end

sol = full(sol);
init_step_sol = sol(1:nx);
other_sol = reshape(sol(nx+1:end), nx+nu+nz, []);

x_unscaled.Fx_opt = other_sol(1,:).'.*Fx_scale;
x_unscaled.rp_opt = other_sol(2,:).'.*rp_scale;
x_unscaled.Fz_opt = other_sol(3,:).'.*Fz_scale;
x_unscaled.Fy_opt = other_sol(4,:).'.*Fy_scale;
x_unscaled.Pow_opt = other_sol(5,:).'.*P_scale;
x_unscaled.u_opt = [init_step_sol(1); other_sol(6,:).'].*u_scale;
x_unscaled.v_opt = [init_step_sol(2); other_sol(7,:).'].*v_scale;
x_unscaled.r_opt = [init_step_sol(3); other_sol(8,:).'].*r_scale;
x_unscaled.ep_opt = [init_step_sol(4); other_sol(9,:).'].*ep_scale;
x_unscaled.ef_opt = [init_step_sol(5); other_sol(10,:).'].*ef_scale;
x_unscaled.alfa_opt = [init_step_sol(6); other_sol(11,:).'].*alfa_scale;

