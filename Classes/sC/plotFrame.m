function plotHandles = plotFrame(M, options)
%
% plots the frame of which position and orientation is defined by the
% homogeneous transform matrix M;
% 
% options: ... 
%
arguments
   M
   options.scale = 1;       % scale for the frame (default to 100% of the unit measure)
   options.parent = gca;    % default parent is the current axes
   options.linewidth = 1.3; % linewidth of the frames' axes
   options.label = [];      % label as pedices of the x, y, z axis names (e.g. <'label', g> gives x_g, y_g, z_g )
   options.fontSize = 12;   % font size of the axes labels
   options.color = [];
end

% extract name valued arguments
scale = options.scale;
parent = options.parent;
linewidth = options.linewidth;
label = options.label;
fontSize = options.fontSize;
color = options.color;

% build the i, j, k versors
x_vec = [M(1:3, 4), M(1:3, 4) + M(1:3, 1).*scale];
y_vec = [M(1:3, 4), M(1:3, 4) + M(1:3, 2).*scale];
z_vec = [M(1:3, 4), M(1:3, 4) + M(1:3, 3).*scale];

frame = line(parent, [x_vec(1, :); y_vec(1, :); z_vec(1, :)].',...
                     [x_vec(2, :); y_vec(2, :); z_vec(2, :)].',...
                     [x_vec(3, :); y_vec(3, :); z_vec(3, :)].',...
                     'linewidth', linewidth);
if isempty(color)                 
    set(frame, {'color'}, {'r';'g';'b'}); % color the axes with r g b order {[1 0 0]; [0 1 0]; [0 0 1]}
else
    set(frame, {'color'}, {color;color;color})
end

% apply the axis labels if the input is given
if ~isempty(label)
    
   labels = text(parent, [x_vec(1,2), y_vec(1,2), z_vec(1,2)],...
                         [x_vec(2,2), y_vec(2,2), z_vec(2,2)],...
                         [x_vec(3,2), y_vec(3,2), z_vec(3,2)],...
                         {strcat('$x_{', label, '}$'), strcat('$y_{', label, '}$'), strcat('$z_{', label, '}$')},...
                         'fontsize', fontSize);

    
end

if nargout == 1
    plotHandles{1} = frame;
    plotHandles{2} = labels;
end
end