function J = BodyJac_localPOE(Goffset, X, q, options)
%
%
%
arguments
    Goffset
    X
    q
    options.EEoffset = eye(4);
end
EEoffset = options.EEoffset;
n = length(q);
J = zeros(6, n, class(q(1)));
Binv = rigidInverse(EEoffset);
for i = n:-1:1
    Binv = Binv*expTw(-X(:, i), q(i));
    J(:,i) = adjoint(Binv)*X(:, i);
    Binv = Binv*rigidInverse(Goffset{i});
end