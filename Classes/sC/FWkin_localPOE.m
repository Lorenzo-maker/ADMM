function [Gglobal, Glocals] = FWkin_localPOE(Goffset, X, q, options)
%
% inputs: Goffset: 1xn cell of 4x4 offset matrices 
%
%         X: 6xn matrix where each row is the unitary twist of each joint
%
%         q: nx1 array where each element is the joint value
%
% options: 'jstart' & 'jend': starting and ending joint for the forward
%                             kinematic computation
%
% outputs the forward kinematics with local POE formulation and each local
% transforms

arguments
    Goffset
    X
    q
    options.jstart = 1
    options.jend = length(q)
end

jstart = options.jstart;
jend = options.jend;

Glocals = cell(1, jend);
Gglobal = eye(4);

for i = jstart:jend
    Glocals{i} = Goffset{i}*expTw(X(:, i), q(i));
    Gglobal = Gglobal*Glocals{i};
end


end