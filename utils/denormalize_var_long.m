function [states_sol, control_sol, algebraic_sol] = denormalize_var_long(X, Nproblems, nx, nu, nz, overlap, dalfa)
%
% x_unscaled = denormalize_var(x_scaled)
%
car_parameters_ocp;

% concatenate solutions
sol = [];
id_consensus = nan(Nproblems-1, 1);
if Nproblems > 0
    index = length(X{1}) - (nx + floor(overlap/2 + 1)*(nx + nu + nz))+(nx+nu+nz+nx);
    sol = cat(1, sol, X{1}(1:index));
    for i = 2:Nproblems-1
        index = length(X{i}) - (nx + floor(overlap/2 + 1)*(nx + nu + nz))+(nx+nu+nz+nx);
        sol = cat(1, sol, X{i}(1 +floor(overlap/2)*(nx+nu+nz)+(nx+nu+nz+nx): index));
    end
    sol = cat(1, sol, X{Nproblems}(1 +floor(overlap/2)*(nx+nu+nz)+(nx+nu+nz+nx): end));
else
    sol = cat(1, sol, X{1});
end
sol = full(sol);
init_step_sol = sol(1:nx);
other_sol = reshape(sol(nx+1:end), nx+nu+nz, []);

states_sol = [init_step_sol, other_sol(nu+nz+1:end, :)];
control_sol = other_sol(1:nu, :);
algebraic_sol = other_sol(nu+1:nu+nz, :);