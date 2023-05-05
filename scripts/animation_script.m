%% animation script


pista = track2casadiSpline(pista); % make casadi splines for faster computation (and animation)
mainAx = pista.full_plot(linspace(0, alfa_end, 5000), 'closed', true);
plot3(traj(1, :), traj(2, :), traj(3, :), 'color', 'y', 'linewidth', 1.5, 'parent', mainAx)
Fig = mainAx.Parent;
Fig.Color = [0.5 0.5 0.5];
Fig.Position = [0,0, 1366, 780];
mainAx.Units = 'normalized';
mainAx.Position = [0, 0, 0.7, 0.7];
hold(mainAx, 'on')
%grid(mainAx, 'off')
mainAx.Visible = 'off';


%steering wheel plot
steering_sol = U_sol(end, :)*10; % assumng 10 steering ratio

steering_axes = axes('parent', Fig, 'Units', 'normalized', 'position', [0.75, 0.15, 0.2, 0.2], 'visible', 'off');
hold(steering_axes, 'on')
%grid(steering_axes, 'off')
axis(steering_axes, 'equal')
steering_transform = hgtransform(steering_axes);
Xsteer = load('Data/coordinateX_steering.mat');
Ysteer = load('Data/coordinateY_steering.mat');
X_patch = load('Data/coordinateX_steering_patch.mat');
Y_patch = load('Data/coordinateY_steering_patch.mat');
Xsteer = Xsteer.X;
Ysteer = Ysteer.Y;
X_patch = X_patch.X_patch;
Y_patch = Y_patch.Y_patch;
surface(Xsteer, Ysteer, Ysteer.*0, 'parent', steering_transform, 'facecolor', 'k');
patch('XData', X_patch, 'YData',Y_patch, 'parent', steering_transform, 'facecolor', 'k');
steering_transform.Matrix = TrotZ(steering_sol(1));

% gauge plot
Gage = uigage(...
    'Parent',Fig, ...
    'Position',[25 550 200  200],...
    'BackColor',[0.35 0.36 0.40],...
    'BackEdgeColor',[0.45 0.46 0.50],...
    'DisplayVisible','off',...
    'DialAngle',[0 90],...
    'DialValues', 0:20:200,...
    'DialVisible','off',...
    'DialLineWidth',2,...
    'DialColor',[0.98 1 1],...
    'TickColor',[0.98 1 1],...
    'TickLineWidth',2,...
    'LabelFontName','Arial',...
    'LabelFontSize',11,...
    'LabelColor',[0.98 1 1],...
    'BandValues',[0 8.2 8.21 8.8 9],...
    'BandColor',[0.80 0.81 0.82 ;...
                 0.80 0.81 0.82 ; ...
                 1.00 0.50    0 ; ...
                 1.00    0    0 ; ...
                 1.00    0    0 ],...
    'UnitString','km/h',...
    'UnitVisible','on',...
    'UnitFontName','Arial',...
    'UnitFontSize',16,...
    'UnitPosition',[0.6 50],...
    'UnitColor',[0.80 0.81 0.82 ],...
    'NeedleLength',1.05);

set(Gage,'NeedlePosition',twist(1, 1)*3.6);

% plot car
car.FWKin_G(X_sol(1,1), X_sol(2:6, 1))
car.plot('plot_track', false, 'parent', mainAx, 'labels', false)

for kk = 1:length(X_sol)-1
   tic
   set(Gage,'NeedlePosition',twist(1, kk)*3.6);
   steering_transform.Matrix = TrotZ(steering_sol(kk));
   time = (1./X_sol(7,kk)).*dalfa;
   % interpolate solutions for smoother animation
   car.FWKin_G(X_sol(1,kk), X_sol(2:6, kk));
   car.updateplot('camera', 1, 'boxSize', [30, 30, 10]);
   drawnow
   t = toc;
   Nframes = floor(time/t); % approx. number of frames avialable for each collocation interval
   
   deltaX = (X_sol(:,kk+1) - X_sol(:,kk))/(Nframes+1);
   delta_steer = (steering_sol(kk+1) - steering_sol(kk))/(Nframes+1);
   for jj =1:Nframes
       steering_transform.Matrix = TrotZ(steering_sol(kk) + delta_steer.*jj);
       car.FWKin_G(X_sol(1,kk) + deltaX(1).*jj, X_sol(2:6,kk) + deltaX(2:6).*jj);
       car.updateplot('camera', 1, 'boxSize', [30, 30, 10]);
       drawnow
   end
   
end