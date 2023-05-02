function [xc_unscaled] = denormalize_xc_long(Xc, d, Nproblems, nx, overlap, dalfa)
%
% x_unscaled = denormalize_var(x_scaled)
%
car_parameters_ocp;

% concatenate solutions
sol = [];
if Nproblems > 0
    sol = cat(1, sol, Xc{1}(1:end - (floor(overlap/2 + 1)*(nx*d))+(nx*d)));
    for i = 2:Nproblems-1

        sol = cat(1, sol, Xc{i}(1 +floor(overlap/2)*(nx*d)+(nx*d): end - (floor(overlap/2 + 1)*(nx*d))+(nx*d)));

    end
    sol = cat(1, sol, Xc{Nproblems}(1 +floor(overlap/2)*(nx*d)+(nx*d): end));
else
    sol = cat(1, sol, Xc{1});
end
sol = full(sol);
xc_unscaled.alfa = sol(1:nx:end)*alfa_scale;
xc_unscaled.ep = sol(2:nx:end)*ep_scale;
xc_unscaled.ef = sol(3:nx:end)*ef_scale;




