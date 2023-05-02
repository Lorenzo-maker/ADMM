function out = subGraph_DersBasisFuns(u, p, U, derOrder, deltaU)
% u = value
% i = knot interval
% p = degree
% U = knot vector


% compute the knot interval (knot span)

i = find(u > U, 1, 'last');

i_prev = i-1;
i_next = i+1;
if u >= U(end)   % check if u value is equal to the last knot
    i = length(U) - p - 1; % the interval is n-p
elseif u <= U(1) % check if u value is equal to the initial knots
    i = p + 1;            % the interval is p + 1
end

if i == p + 1
    i_prev = p + 1;
    i_next = i+1;
end

if i == length(U) - p - 1
    i_next = length(U) - p - 1;
    i_prev = i-1;
end

U_middle(1) = U(i);
U_middle(2) = U(i + 1);
u_num = u;
import casadi.*
u = SX.sym('u');

N_middle = derBasisFun_i(u, i, p, U, derOrder);
N_prev = derBasisFun_i(u, i_prev, p, U, derOrder);
N_next = derBasisFun_i(u, i_next, p, U, derOrder);

%% computing symbolic if_else

% reduce the computational graph as low as possible, yet mantaining
% numerical stability by considering neighbour knot spans when delta alpha
% might fall in there during optimization
if u_num - deltaU <= U_middle(1)
    out = if_else(u < U_middle(1), N_prev, N_middle);
elseif u_num + deltaU >= U_middle(2)
    out = if_else(u < U_middle(2), N_middle, N_next);
else
    out = N_middle;
end
if u_num - deltaU < U_middle(1) && u_num + deltaU >= U_middle(2)
    out = if_else(u >= U_middle(1) & u <= U_middle(2), N_middle,...
        if_else(u < U_middle(1), N_prev, N_next ));
end

out = Function('N', {u}, {out'});

end