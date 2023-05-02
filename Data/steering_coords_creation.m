%% 
% Define the steering wheel center and radius
center = [0 0];
radius = 1;

% Define the number of spokes for the wheel
spokes = 50;

% Define the angle for each spoke
theta = (2*pi)/spokes;

% Pre-allocate arrays to store the x and y coordinates of each spoke
x_spokes = zeros(spokes, 1);
y_spokes = zeros(spokes, 1);

x_spokes_internal = zeros(spokes, 1);
y_spokes_internal = zeros(spokes, 1);
% Generate the x and y coordinates of each spoke
for i = 1:spokes+1
    x_spokes(i) = center(1) + radius*cos((i-1)*theta);
    y_spokes(i) = center(2) + radius*sin((i-1)*theta);
    
    x_spokes_internal(i) = center(1) + (radius- 0.2)*cos((i-1)*theta);
    y_spokes_internal(i) = center(2) + (radius- 0.2)*sin((i-1)*theta);
end
X = [x_spokes, x_spokes_internal];
Y = [y_spokes, y_spokes_internal];
X_patch = [-0.8 0.8 0.8 0.34 0.15 0.15 -0.15 -0.15 -0.34 -0.8 -0.8];
Y_patch = [0.1 0.1 -0.1 -0.1 -0.3 -0.8 -0.8 -0.3 -0.1 -0.1 0.1];
save('coordinateX_steering.mat', 'X')
save('coordinateY_steering.mat', 'Y')
save('coordinateX_steering_patch.mat', 'X_patch')
save('coordinateY_steering_patch.mat', 'Y_patch')
% Create a patch object to represent the steering wheel
patch_handle = surface([x_spokes, x_spokes_internal], [y_spokes, y_spokes_internal], [y_spokes, y_spokes_internal].*0);
patch(X_patch, Y_patch, 'k')
% Set the patch properties to create a 3D-like appearance
set(patch_handle, 'EdgeColor', 'k', 'LineWidth', 2, 'FaceColor', 'k');

% Limit the axis and add a grid to the plot
axis equal
grid on