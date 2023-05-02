function animateRotation(V, F, integrationSol, exactSol, options)
%
% animateRotation(V, F, integrationSol, exactSol, options)
%
% function that builds the graphical objects tho animate a rotated body
%   
% V: vertex of the body to rotate
%
% F: faces of the body to rotate
%
% integrationSol: 3x3xN matrix representing the integrated rotation matrix
%       at each discretization time interval
%
% exactSol: 3x3xN matrxi representing the exact solution of the  integrated
%       rotation matrix at each discretization time interval
% 
% options: alphaExact -> transparency parameter for the body rotated with
%                        the exact solution (number in [0,1])
%
%          alphaIntegrated -> transparency parameter for the body rotated
%                           with the numerical integration solution
%                           (number in [0,1])

arguments
    V
    F
    integrationSol
    exactSol
    options.alphaExact = 0.5;
    options.alphaIntegrated = 0.5;
end

alphaExact = options.alphaExact;
alphaIntegrated = options.alphaIntegrated;

% initialize figure
figure; hold on; axis equal; box on; grid on; warning off;

transform = hgtransform(gca);            % transform handle for exact rotation
transform_integrated = hgtransform(gca); % transform handle for numerical integrated rotation

% global frame drawing
plotFrame(eye(4), 'parent', gca, 'scale', 3, 'linewidth', 3, 'label', 'g', 'fontSize', 12);

% local frame drawing
local_frame = plotFrame(eye(4), 'parent', transform, 'scale', 1.5, 'label', 'l', 'fontSize', 12);

% object drawing with exact solution
obj_exact = patch(transform, 'Vertices', V, 'Faces', F, 'FaceColor',"b", "EdgeColor",'k', 'facealpha', alphaExact, 'edgealpha', alphaExact, 'linewidth', 1.1);

% object drawing with numerical integration
obj_num = patch(transform_integrated, 'Vertices', V, 'Faces', F, 'FaceColor',"green", "EdgeColor",'k', 'facealpha', alphaIntegrated, 'edgealpha', alphaIntegrated);

% some axis options
set(gca, 'xlim', [-3 3], 'ylim', [-3 3], 'zlim', [-3 3]);
xlabel('x')
ylabel('y')
zlabel('z')
view([45 45])

%% animation method 1 using the hgtransform handles
% for i = 1:size(exactSol, 3)
%     
%     % animation method 1 by applying homogeneous matrix to transform handle
%     transform_integrated.Matrix = RpTohomogeneous(integrationSol(:, :, i), [0; 0; 0]);
%     transform.Matrix = RpTohomogeneous(exactSol(:, :, i), [0; 0; 0]);
% 
%     pause(0.01)
% end

%% animation method 2 by transforming directly the vertices

Vexact = obj_exact.Vertices.';
Vnum = obj_num.Vertices.';
[poslabelX, poslabelY, poslabelZ] = local_frame{2}.Position;
posAx = [[local_frame{1}.XData]; [local_frame{1}.YData]; [local_frame{1}.ZData]];

for i = 1:size(exactSol, 3)

    exactMat = exactSol(:, :, i);
    axpt = exactMat*posAx;                    % compute updated axes point

    local_frame{1}(1).XData = axpt(1,1:2);    % update axes points
    local_frame{1}(2).XData = axpt(1,3:4);
    local_frame{1}(3).XData = axpt(1,5:6);
    local_frame{1}(1).YData = axpt(2,1:2);
    local_frame{1}(2).YData = axpt(2,3:4);
    local_frame{1}(3).YData = axpt(2,5:6);
    local_frame{1}(1).ZData = axpt(3,1:2);
    local_frame{1}(2).ZData = axpt(3,3:4);
    local_frame{1}(3).ZData = axpt(3,5:6);

    labelpos{1} = (exactMat*poslabelX.').';   % compute updated label positions
    labelpos{2} = (exactMat*poslabelY.').';
    labelpos{3} = (exactMat*poslabelZ.').';

    [local_frame{2}.Position] = labelpos{:};  % update label position

    obj_exact.Vertices = (exactMat*Vexact).'; % update object vertices with exact results

    obj_num.Vertices = (integrationSol(:, :, i)*Vnum).'; % update object vertices with numerical integration results

    pause(0.01)  % update graphics
end

end