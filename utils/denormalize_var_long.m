function [states_sol, control_sol, algebraic_sol] = denormalize_var_long(X, Nproblems, nx, nu, nz, id)
%
% x_unscaled = denormalize_var(x_scaled)
%
%car_parameters_ocp;

% concatenate solutions
sol = [];


id_consensus = nan(Nproblems-1, 1);
if Nproblems > 0
    sol = cat(1, sol, X{1}(1:(id.t{1}(1)-1)*(nx+nu+nz)));
    for i = 2:Nproblems-1                
        sol = cat(1, sol, X{i}((id.h{i}(1)-1)*(nx+nu+nz)+1:(id.t{i}(1)-1)*(nx+nu+nz)));
    end
    sol = cat(1, sol, X{Nproblems}((id.h{Nproblems}(1)-1)*(nx+nu+nz)+1: end));
else
    sol = cat(1, sol, X{1});
end

sol = full(sol);
init_step_sol = sol(1:nx);
other_sol = reshape(sol(nx+1:end), nx+nu+nz, []);

states_sol = [init_step_sol, other_sol(nu+nz+1:end, :)];
control_sol = other_sol(1:nu, :);
algebraic_sol = other_sol(nu+1:nu+nz, :);


