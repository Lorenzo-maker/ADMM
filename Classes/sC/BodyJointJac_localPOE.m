function J = BodyJointJac_localPOE(Goffset, X, q, j)
%
% returns the body jacobian w.r.t. to the selected joint j.
% basically the returned jacobian maps to twists of the kinematic chain
% expressed in the j-th components and origin (polo)
%

n = length(q);
J = zeros(6, n, class(q(1)));

% compute first the j-th joint jacobian contribution
B = eye(4);

for i = j:-1:1 % compute backwards the kinematic chain
   B = B*expTw(-X(:, i), q(i));
   J(:, i) = adjoint(B)*X(:, i);
   B = B*rigidInverse(Goffset{i});
end

B = expTw(X(:, j), q(j))*Goffset{j};

for i = j+1:n % compute forwards the kinematic chain
    B = B*expTw(X(:, i), q(i));
    J(:,i) = adjoint(B)*X(:, i);
    B = B*Goffset{i};
end