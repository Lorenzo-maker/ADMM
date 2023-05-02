%% UIGAGE user interface needle gage
%
% The class *uigage* can be used in a GUI to display outputs that are suitable 
% for a needle gage type display. It behaves in a very similar fashion to 
% other Matlab *uicontrol* objects.
% Many configuration properties are available which can be accessed through
% SET/GET methods. The examples below show two very different gages that
% can be generated using *uigage*.


%% Example 1: Tachometer Gage
% Creates an engine tacho gage assigning properties during the object
% creation and moving the gage needle to the desired value.
% (Using *SET/GET* methods to modify and retrieve properties) 

% Creating figure and modifying some of its properties
Fig1 = figure;
FigPos = get(Fig1,'Position');
set(Fig1,...
    'Position',[FigPos(1:2) 600 400],...
    'Name','Tacho Example',...
    'Color',[0.2 0.21 0.22]);

% Creating tacho gage
TachoGage = uigage(...
    'Parent',Fig1, ...
    'Position',[100 1 400 400],...
    'BackColor',[0.35 0.36 0.40],...
    'BackEdgeColor',[0.45 0.46 0.50],...
    'DisplayVisible','off',...
    'DialAngle',[0 90],...
    'DialValues',0:9,...
    'DialVisible','off',...
    'DialLineWidth',2,...
    'DialColor',[0.98 1 1],...
    'TickColor',[0.98 1 1],...
    'TickLineWidth',2,...
    'LabelFontName','Arial',...
    'LabelFontSize',16,...
    'LabelColor',[0.98 1 1],...
    'BandValues',[0 8.2 8.21 8.8 9],...
    'BandColor',[0.80 0.81 0.82 ;...
                 0.80 0.81 0.82 ; ...
                 1.00 0.50    0 ; ...
                 1.00    0    0 ; ...
                 1.00    0    0 ],...
    'UnitString','x1000 rpm',...
    'UnitVisible','on',...
    'UnitFontName','Arial',...
    'UnitFontSize',16,...
    'UnitPosition',[0.6 50],...
    'UnitColor',[0.80 0.81 0.82 ],...
    'NeedleLength',1.05);

% Changing needle position
set(TachoGage,'NeedlePosition',7.75);


%% Example 2: Pressure Gage
% Creates a pressure gage assigning properties during the object
% creation and moving the gage needle to the desired value.
% Additional properties are modified after the object is created.
% (Using *OBEJCT.PROPERTY* approach to modify and retrieve properties) 

% Creating figure and modifying some of its properties
Fig2 = figure;
FigPos = Fig2.Position;
Fig2.Position = [FigPos(1) FigPos(2)-240 400 200];
Fig2.Name = 'Pressure Gage Example';
Fig2.NumberTitle = 'off';
Fig2.MenuBar = 'none';
Fig2.Color = [0.97 0.99 0.99];

% Creating pressure gage
PressureGage = uigage(Fig2, ...
    'Position',[50 -125 300 300],...
    'DialValues',1:5,...
    'DialAngle',[-120 120],...
    'DisplayVisible','off',...
    'BackColor','none',...
    'BackEdgeColor','none',...
    'LabelDir','out',...
    'LabelFontName','Arial',...
    'LabelFontSize',20,...
    'LabelFontWeight','normal',...
    'BandPosition',[0.8 1],...
    'BandValues',[1 3.59 3.6 3.99 4 5],...
    'BandColor',[0.5725 0.8157 0.3137 ;...
                 0.5725 0.8157 0.3137 ;...
                 1.0000 1.0000      0 ;...
                 1.0000 1.0000      0 ;...
                 1.0000      0      0 ;...
                 1.0000      0      0 ],...
    'TickPosition',0.9,...
    'TickDir','center',...
    'TickLength',0.2,...
    'MinorTickVisible','off',...
    'KnobSize',0.12,...
    'KnobColor',[0.1 0.11 0.11],...
    'NeedleLength',1.0,...
    'NeedleColor',[0.1 0.11 0.11],...
    'NeedleWidth',[0.03 0.01]);

% Changing needle position
PressureGage.NeedlePosition = 3.85;

% Modifying additional properties
PressureGage.UnitVisible = 'on';
PressureGage.UnitString = 'bar';
PressureGage.UnitFontName = 'Arial';
PressureGage.UnitFontSize = 20;



