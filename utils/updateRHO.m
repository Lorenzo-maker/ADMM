function newRHO = updateRHO(oldRHO, X, Z_current, Z_previous, Y, ADMM_iteration)
%
%
%

tauMax = 1.5;

%primal_res = vecnorm((X - Z_current));
primal_res = vecnorm((X - Z_previous));
dual_res = vecnorm((Z_current - Z_previous));

ratio = sqrt( vecnorm((X - Z_current))./ vecnorm(oldRHO.*(Z_current - Z_previous)) );
ratio2 = sqrt( vecnorm(oldRHO.*(Z_current - Z_previous)) ./ vecnorm((X - Z_current)) );

if ratio < tauMax && ratio >= 1
    tau = ratio;
elseif ratio < 1 && ratio > 1/tauMax
    tau = ratio2;
else
    tau = tauMax;
end

tau = 2;
mu = 10;
if (primal_res > dual_res*mu) && (ADMM_iteration > 0)
    newRHO = oldRHO.*tau;
    newRHO = min(newRHO, 1000);
    return
end

if (mu*primal_res < dual_res) && (ADMM_iteration > 0)
    newRHO = oldRHO./tau;
    newRHO = max(newRHO, 1);
    return
end

newRHO = 1.01*oldRHO;
end