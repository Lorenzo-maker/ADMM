classdef uigage < matlab.mixin.SetGet
%  UIGAGE Create user interface needle gage.
%     H = uigage('PropertyName1',value1,'PropertyName2',value2,...) 
%     creates a needle gage user interface control in the current figure 
%     window and returns a handle to it in H.
%  
%     H = uigage(FIG,...) creates a needle gage user interface control in the
%     specified figure.
%  
%     uigage properties can be set at object creation time using
%     PropertyName/PropertyValue pair arguments to uigage, or 
%     changed later using the SET command.
%  
%     Execute GET(H) to see a list of uigage object properties and
%     their current values. Execute SET(H) to see a list of uigage
%     object properties and legal property values.
%  
%     Examples:
%        Example 1:
%             %creates default uigage in a new figure
%             H = uigage; 
%       
%        Example 2:
%             %creates a figure with a uigage in it and moves the needle
%             %some properties are modified after the object is created
%             Fig1 = figure;
%             H = uigage(Fig1); 
%             set(H,'NeedlePosition',0.7)
%             set(H,...
%                 'BandValues',[0.9 0.95 1.0],...
%                 'BandPosition',[0.95 1],...
%                 'DialMinorValues',0:0.05:1,...
%                 'TickLength',0.2,...
%                 'TickPosition',1.05,...
%                 'LabelDir','out');
%
%  See also SET, GET
  
%
% Class developed by Eduardo Nigro in September of 2017
% Version 1.00
%
% Updates:
%
% Version 1.01
% - Sep 2017 - Fixed 'Visible' property bug
% - Sep 2017 - Made Parent property immutable
% - Sep 2017 - Added Enable property
%
% Version 1.02
% - Nov 2017 - Fixed bug for gage creation with no input arguments
% - Nov 2017 - Fixed bug for very small separation between color band values
% - Nov 2017 - Fixed bug with property assignment sequence in constructor
%              to handle special cases of LabelString property
%
% Version 1.03
% - Mar 2019 - Added uipanel as a possible parent
%            - Fixed default values of OffBandColor (to match default band size)
%
% Version 1.04
% - Sep 2020 - Using tranformation group for faster needle position update
%            - Improved speed of gage numerical display update
%

    % Internal properties (angle arrays)
    properties (Access = private, Hidden = true)
        
        ThtDial                     % Dial ring
        ThtBand                     % Dial color band
        ThtTick                     % Dial major tick marks
        ThtMinorTick                % Dial minor tick marks
        
        % Disabled component colors
        OffBackColor              = [0.9 0.9 0.9];
        OffBackEdgeColor          = [0.75 0.75 0.75];
        OffBandColor              = 0.85*ones(3,3);
        OffDialColor              = [0.6 0.6 0.6];
        OffLabelColor             = [0.65 0.65 0.65];
        OffTickColor              = [0.6 0.6 0.6];
        OffMinorTickColor         = [0.65 0.65 0.65];
        OffKnobColor              = [0.65 0.65 0.65];
        OffKnobEdgeColor          = [0.6 0.6 0.6];
        OffNeedleColor            = [0.65 0.65 0.65];
        OffNeedleEdgeColor        = [0.6 0.6 0.6];
        OffDisplayForegroundColor = [0.6 0.6 0.6];
        OffDisplayBackgroundColor = [0.7 0.7 0.7];
        OffDisplayEdgeColor       = [0.5 0.5 0.5];
        OffUnitColor              = [0.6 0.6 0.6];
        
    end
    
    % Gage components (handles)
    properties (Access = private, Hidden = true)
        Axis                        % Axis used to build gage
        Back                        % Back panel
        Dial                        % Dial ring
        Band                        % Color band
        Label                       % Dial labels
        Tick                        % Dial major tick marks
        MinorTick                   % Dial major tick marks
        Knob                        % Needle knob
        Needle                      % Needle
        NeedleGrp                   % Needle transf. group
        Display                     % Numerical display
        Unit                        % Units display
    end
    
    properties (SetAccess = immutable)
        Parent                      % Parent figure handle
    end
    
    % Dependent properties
    properties (Dependent)
        
        Box                         % uigage delimiting box :  { 'on'  'off' }
        Position                    % uigage position :  [ XLowerLeftCorner  YLowerLeftCorner  Width  Height ]
        Units                       % uigage units :  { 'pixels'  'normalized' }
        Visible                     % uigage visibility :  { 'on'  'off' }
        Enable                      % uigage enable : {'on'  'off'}

        BackColor                   % back plane face color :  [ R  G  B ]  triplet
        BackEdgeColor               % back plane edge color :  [ R  G  B ]  triplet
        BackLineWidth               % back plane edge line width
        BackVisible                 % back plane visibility :  { 'on'  'off' }
        
        BandValues                  % color band values :  [ Value1  Value2  ...  ValueN ] (1xN)
        BandColor                   % color band colors :  [ R1  G1  B1 ;   R1  G1  B1 ;  ...   RN  GN  BN ] (Nx3)
        BandPosition                % color band position :  [ Radius1  Radius2 ]
        BandVisible                 % color band visibility :  { 'on'  'off' }
        
        DialAngle                   % dial angle range (ZERO correponds to -Y axis) :  [ Angle1  Angle2 ] 
        DialValues                  % dial major tick values :  [ Value1  Value2  ...  ValueN ] (1xN)
        DialMinorValues             % dial minor tick values :  [ Value1  Value2  ...  ValueM ] (1xM)
        DialColor                   % dial face color :  [ R  G  B ]  triplet
        DialLineWidth               % dial line width
        DialVisible                 % dial visibility :  { 'on'  'off' }

        LabelString                 % label string :  {'Name1' ,  'Name2' ,  ... , 'NameN'}
        LabelFontName               % label font name
        LabelFontSize               % label font size
        LabelFontWeight             % label font weight :  { 'normal'  'bold' }
        LabelColor                  % label color :  [ R  G  B ]  triplet
        LabelDir                    % label position with respect to dial :  { 'in'  'center'  'out' }
        LabelOrientation            % label orientation :  {'horizontal'  'radial'  'tangential'}
        LabelVisible                % label visibility :  { 'on'  'off' }
        
        TickColor                   % major tick color :  [ R  G  B ]  triplet
        TickDir                     % major tick position with respect to dial :  { 'in'  'center'  'out' }
        TickLength                  % major tick length with repsect to unit radius 
        TickLineWidth               % major tick line width
        TickPosition                % major tick radial position :  Radius
        TickVisible                 % major tick visibility :  { 'on'  'off' }
        
        MinorTickColor              % minor tick color :  [ R  G  B ]  triplet
        MinorTickDir                % minor tick position with respect to dial :  { 'in'  'center'  'out' }
        MinorTickLength             % minor tick length with repsect to unit radius 
        MinorTickLineWidth          % minor tick line width
        MinorTickPosition           % minor tick radial position :  Radius
        MinorTickVisible            % minor tick visibility :  { 'on'  'off' }
        
        KnobSize                    % knob diameter with repsect to unit radius 
        KnobColor                   % knob color :  [ R  G  B ]  triplet
        KnobEdgeColor               % knob edge color :  [ R  G  B ]  triplet
        KnobLineWidth               % knob edge line width
        
        NeedleColor                 % needle color :  [ R  G  B ]  triplet
        NeedleEdgeColor             % needle edge color :  [ R  G  B ]  triplet
        NeedleLineWidth             % needle edge line width :  [ ValueBase ]  or  [ ValueBase  ValueTip ]
        NeedleLength                % needle length with repsect to unit radius 
        NeedleWidth                 % needle width with repsect to unit radius 
        NeedleOffset                % needle offset on knob with repsect to unit radius 
        % needle shape points :  [ X1  Y1 ;  X2  Y2 ;  ...  XN  YN ] (Nx2)
        % - superseeds 'NeedleLength' and 'NeedleWidth' if defined by user
        % - automatically calculated if 'NeedleLength' or 'NeedleWidth' are set
        NeedlePoints    
        NeedlePosition              % needle position :  Value
        
        DisplayFontName             % numeric display font name
        DisplayFontSize             % numeric display font size
        DisplayFontWeight           % numeric display font weight :  { 'normal'  'bold' }
        DisplayForegroundColor      % numeric display foreground color :  [ R  G  B ]  triplet
        DisplayBackgroundColor      % numeric display background color :  [ R  G  B ]  triplet
        DisplayEdgeColor            % numeric display edge color :  [ R  G  B ]  triplet
        DisplayLineWidth            % numeric display edge line width
        DisplayPosition             % numeric display position :  [ Radius  Angle  BoxWidth  BoxHeight ]
        DisplayVisible              % numeric display visibility :  { 'on'  'off' }
        
        UnitString                  % unit string
        UnitFontName                % unit font name
        UnitFontSize                % unit font size
        UnitFontWeight              % unit font weight :  { 'normal'  'bold' }
        UnitColor                   % unit font color :  [ R  G  B ]  triplet
        UnitPosition                % unit position :  [ Radius  Angle ]
        UnitVisible                 % unit visibility :  { 'on'  'off' }
        
    end
    
    % Corresponding private properties used for SET/GET methods
    properties (Access = private)
        
        PrvtBox = 'off'
        PrvtPosition = []
        PrvtUnits = 'pixels'
        PrvtVisible = 'on'
        PrvtEnable = 'on';
        
        PrvtBackColor = [0.95 0.96 0.96]
        PrvtBackEdgeColor = [0.85 0.86 0.86]
        PrvtBackLineWidth = 0.5
        PrvtBackVisible = 'on'
        
        PrvtBandValues = [0.8 0.9 1.0]
        PrvtBandColor = [1 1 0; 1 0.5 0; 1 0 0]
        PrvtBandPosition = [0.9 1.0]
        PrvtBandVisible = 'on'

        PrvtDialAngle = [-60 60]
        PrvtDialValues = 0:0.1:1
        PrvtDialMinorValues = 0:0.02:1
        PrvtDialColor = [0.20 0.21 0.21]
        PrvtDialLineWidth = 1.5
        PrvtDialVisible = 'on'

        PrvtLabelString = ''
        PrvtLabelFontName = 'Helvetica'
        PrvtLabelFontSize = 10
        PrvtLabelFontWeight = 'bold'
        PrvtLabelColor = [0.25 0.26 0.26]
        PrvtLabelDir = 'in'
        PrvtLabelOrientation = 'horizontal'
        PrvtLabelVisible = 'on'
        
        PrvtTickColor = [0.20 0.21 0.21]
        PrvtTickDir = 'in'
        PrvtTickLength = 0.15
        PrvtTickLineWidth = 1.5
        PrvtTickPosition = 1
        PrvtTickVisible = 'on'

        PrvtMinorTickColor = [0.25 0.26 0.26]
        PrvtMinorTickDir = 'in'
        PrvtMinorTickLength = 0.1
        PrvtMinorTickLineWidth = 0.5
        PrvtMinorTickPosition = 1
        PrvtMinorTickVisible = 'on'
        
        PrvtKnobSize = 0.16
        PrvtKnobColor = [0.25 0.26 0.26]
        PrvtKnobEdgeColor = [0.50 0.51 0.51]
        PrvtKnobLineWidth = 0.5
        
        PrvtNeedleColor = [0.9 0.05 0.1]
        PrvtNeedleEdgeColor = [0.50 0.51 0.51]
        PrvtNeedleLineWidth = 0.5
        PrvtNeedleLength = 0.93
        PrvtNeedleWidth = [0.08 0.01]
        PrvtNeedleOffset = -0.1
        PrvtNeedlePoints
        PrvtNeedlePosition = 0
        
        PrvtDisplayFontName = 'Helvetica'
        PrvtDisplayFontSize = 12
        PrvtDisplayFontWeight = 'normal'
        PrvtDisplayForegroundColor = [0.20 0.21 0.21]
        PrvtDisplayBackgroundColor = [0.90 0.91 0.91]
        PrvtDisplayEdgeColor = [0.50 0.51 0.51]
        PrvtDisplayLineWidth = 1
        PrvtDisplayPosition = [0.6 0 0.5 0.2]
        PrvtDisplayVisible = 'on'
        
        PrvtUnitString = ' '
        PrvtUnitFontName = 'Helvetica'
        PrvtUnitFontSize = 10
        PrvtUnitFontWeight = 'normal'
        PrvtUnitColor = [0.20 0.21 0.21]
        PrvtUnitPosition = [0.5 180]
        PrvtUnitVisible = 'off'

    end
    
    methods
        
% =========================================================================
%                        CONSTRUCTOR METHOD
% =========================================================================           

        function obj = uigage(varargin)
            
            % UIGAGE uses several graphic objects to create a needle gage
            % uicontrol. The radial positions of different components are 
            % based on the gage unit radius.
            %
            % As a passive uicontrol, uigage only has properties but no
            % public methods. For simplicity it does not contain any 
            % callback functionality.
            %
            % Use the 'NeedlePosition' property to continuously update 
            % the needle position inside your application.
            
            % Checking possible contructor inputs
            if nargin == 0
                
                % Assigning parent figure
                obj.Parent = gcf;
                
                % Assinging start index
                iStart = 0;
                    
            else
                
                % Checking for odd number of inputs
                if rem(nargin,2) ~= 0
                    
                    % Checking for figure handle
                    if ishandle(varargin{1})
                        HandleAux = varargin{1};
                        if strcmp(HandleAux.Type,'figure') || strcmp(HandleAux.Type,'uipanel')
                            obj.Parent = HandleAux;
                        else
                            error('Handle does not appear to be a figure or panel')
                        end
                    else
                        error('First input argument must be a figure handle')
                    end
                    
                    % Assinging start index
                    iStart = 2;
                    
                else
                    
                    % Assinging start index
                    iStart = 1;
                    
                end
                
                % Looping through property (short list)
                for i = iStart:2:nargin-1
                    switch lower(varargin{i})
                        case 'parent'
                            obj.Parent = varargin{i+1};
                        case 'position'
                            obj.Position = varargin{i+1};
                        case 'units'
                            obj.Units = varargin{i+1};
                    end
                end
                
            end
            
            % Assigning parent if needed
            if isempty(obj.Parent)
                obj.Parent = gcf;
            end
            
            % Assigning position if needed
            if isempty(obj.Position)
                FigPos = obj.Parent.Position;
                GageSize = 0.8*min(FigPos(3:4));
                obj.Position = [...
                    0.5*(FigPos(3)-GageSize) ...
                    0.5*(FigPos(4)-GageSize) ...
                    GageSize ...
                    GageSize];
            end
            
            % Assigning units if needed
            if isempty(obj.Units)
                obj.Units = obj.Parent.Units;
            end
            
            % Creating base axis for gage
            obj.Axis = axes(...
                'Parent',obj.Parent,...
                'Units',obj.Units,...
                'Position',obj.Position,...
                'Box','on',...
                'NextPlot','add',...
                'Color','none',...
                'XLim',[-1.2 1.2],...
                'YLim',[-1.2 1.2],...
                'XColor','none',...
                'YColor','none',...
                'XTickLabel',[],...
                'YTickLabel',[],...
                'DeleteFcn',{@obj.deleteaxis});
            
            % Calculating angle arrays used in by gage components
            [obj.ThtDial,obj.ThtBand,obj.ThtTick,obj.ThtMinorTick] = obj.calculateAngles(...
                obj.DialAngle,obj.BandValues,obj.DialValues,obj.DialMinorValues);
            
            % Creating gage components from bottom to top
            obj.drawBack
            obj.drawBand
            obj.drawDial
            obj.drawMinorTick
            obj.drawTick
            obj.drawLabel
            obj.drawDisplay
            obj.drawUnit
            obj.drawKnob
            obj.drawNeedle
            
            % Assigning user defined properties (long list)
            if nargin ~= 0
                obj.putProperties(nargin,varargin,iStart)
            end
            
            % Looping through property input so special
            % case properties initialize last
            % (Property assignment sequence matters for these folks...)
            for i = iStart:2:nargin-1
                switch lower(varargin{i})
                    case 'labelstring'
                        obj.LabelString = varargin{i+1};
                    case 'enable'
                        obj.Enable = varargin{i+1};
                end
            end

            % Updating needle position
            obj.moveNeedle(obj.NeedlePosition);
            
        end
        
% =========================================================================           
%                       DEFINING SET/GET METHODS    
% =========================================================================           
        
% SET/GET METHODS FOR GAGE MAIN PROPERTIES

        function set.Box(obj,Value)
            
            % Checking input values and updating private property
            if ~strcmpi(Value,'on') &&  ~strcmpi(Value,'off')
                error('''Box'' property must be ''on'' or ''off''');
            end
            obj.PrvtBox = Value;
            
            % Updating relevant gage components
            if strcmpi(obj.PrvtBox,'on') 
                obj.Axis.XColor = [0.50 0.51 0.51];
                obj.Axis.YColor = [0.50 0.51 0.51];
            else
                obj.Axis.XColor = 'none';
                obj.Axis.YColor = 'none';
            end
            
        end
        function Value = get.Box(obj)
            Value = obj.PrvtBox;
        end
        
        function set.Position(obj,Value)

            % Checking input values and updating private property
            if size(Value,1) ~= 1  ||  size(Value,2) ~= 4
                error('Size of position array must be (1,4)');
            end
            obj.PrvtPosition = Value;
        
            % Updating relevant gage components
            obj.Axis.Position = obj.PrvtPosition;
            
        end
        function Value = get.Position(obj)
            Value = obj.PrvtPosition;
        end

        function set.Units(obj,Value)
         
            % Checking input values and updating private property
            if ~strcmpi(Value,'pixels') &&  ~strcmpi(Value,'normalized')
                error('''Units'' property must be ''pixels'' or ''normalized''');
            end
            obj.PrvtUnits = Value;
            
            % Updating relevant gage components
            obj.Axis.Units = obj.PrvtUnits;
            
        end
        function Value = get.Units(obj)
            Value = obj.PrvtUnits;
        end

        function set.Visible(obj,Value)
            
            % Checking input values and updating private property
            if ~strcmpi(Value,'on') &&  ~strcmpi(Value,'off')
                error('''Visible'' property must be ''on'' or ''off''');
            end
            obj.PrvtVisible = Value;
            
            % Updating relevant gage components
            obj.Axis.Visible = obj.PrvtVisible;
            obj.Back.Visible = obj.PrvtVisible;
            obj.Knob.Visible = obj.PrvtVisible;
            obj.Needle.Visible = obj.PrvtVisible;
            obj.Dial.Visible = obj.PrvtVisible;
            obj.Unit.Visible = obj.PrvtVisible;
            for i = 1:length(obj.Display)
                obj.Display{i}.Visible = obj.PrvtVisible;
            end
            for i = 1:length(obj.Band)
                obj.Band{i}.Visible = obj.PrvtVisible;
            end
            for i = 1:length(obj.Label)
                obj.Label{i}.Visible = obj.PrvtVisible;
            end
            for i = 1:length(obj.Tick)
                obj.Tick{i}.Visible = obj.PrvtVisible;
            end
            for i = 1:length(obj.MinorTick)
                obj.MinorTick{i}.Visible = obj.PrvtVisible;
            end
            
        end
        function Value = get.Visible(obj)
            Value = obj.PrvtVisible;
        end
        
        function set.Enable(obj,Value)
            
            % Checking input values updating private property
            if ~strcmpi(Value,'on') &&  ~strcmpi(Value,'off')
                error('''Enable'' property must be ''on'' or ''off''');
            end
            
            % Checking for change in value
            if ~strcmpi(obj.PrvtEnable,Value)
                
                % Updating private property
                obj.PrvtEnable = Value;
                
                % Checking for disable state
                if strcmpi(obj.PrvtEnable,'on')
                    
                    % Applying original colors
                    obj.Back.FaceColor = obj.PrvtBackColor;
                    obj.Back.EdgeColor = obj.PrvtBackEdgeColor;
                    obj.Dial.Color = obj.PrvtDialColor;
                    for i = 1:length(obj.Label)
                        obj.Label{i}.Color = obj.PrvtLabelColor;
                    end
                    for i = 1:length(obj.Tick)
                        obj.Tick{i}.Color = obj.PrvtTickColor;
                    end
                    for i = 1:length(obj.MinorTick)
                        obj.MinorTick{i}.Color = obj.PrvtMinorTickColor;
                    end
                    obj.Knob.FaceColor = obj.PrvtKnobColor;
                    obj.Knob.EdgeColor = obj.PrvtKnobEdgeColor;
                    obj.Needle.FaceColor = obj.PrvtNeedleColor;
                    obj.Needle.EdgeColor = obj.PrvtNeedleEdgeColor;
                    obj.Display{2}.Color = obj.PrvtDisplayForegroundColor;
                    obj.Display{1}.FaceColor = obj.PrvtDisplayBackgroundColor;
                    obj.Display{1}.EdgeColor = obj.PrvtDisplayEdgeColor;
                    obj.Unit.Color = obj.PrvtUnitColor;
                    
                else
                    
                    % Applying grayed out colors
                    obj.Back.FaceColor = obj.OffBackColor;
                    obj.Back.EdgeColor = obj.OffBackEdgeColor;
                    obj.Dial.Color = obj.OffDialColor;
                    for i = 1:length(obj.Label)
                        obj.Label{i}.Color = obj.OffLabelColor;
                    end
                    for i = 1:length(obj.Tick)
                        obj.Tick{i}.Color = obj.OffTickColor;
                    end
                    for i = 1:length(obj.MinorTick)
                        obj.MinorTick{i}.Color = obj.OffMinorTickColor;
                    end
                    obj.Knob.FaceColor = obj.OffKnobColor;
                    obj.Knob.EdgeColor = obj.OffKnobEdgeColor;
                    obj.Needle.FaceColor = obj.OffNeedleColor;
                    obj.Needle.EdgeColor = obj.OffNeedleEdgeColor;
                    obj.Display{2}.Color = obj.OffDisplayForegroundColor;
                    obj.Display{1}.FaceColor = obj.OffDisplayBackgroundColor;
                    obj.Display{1}.EdgeColor = obj.OffDisplayEdgeColor;
                    obj.Unit.Color = obj.OffUnitColor;
                    
                end
                
                % Redrawing band (special case)
                obj.drawBand('updated')
                
            end
            
        end
        function Value = get.Enable(obj)
            Value = obj.PrvtEnable;
        end
        
        
% SET/GET METHODS FOR GAGE BACKGROUND

        function set.BackColor(obj,Value)
           
            % Checking input values and updating private property
            obj.PrvtBackColor = Value;

            % Updating relevant gage components
            if strcmpi(obj.Enable,'on')
                obj.Back.FaceColor = obj.PrvtBackColor;
            else
                obj.Back.FaceColor = obj.OffBackColor;
            end
            
        end
        function Value = get.BackColor(obj)
            Value = obj.PrvtBackColor;
        end
        
        function set.BackEdgeColor(obj,Value)

            % Checking input values and updating private property
            obj.PrvtBackEdgeColor = Value;
            
            % Updating relevant gage components
            if strcmpi(obj.Enable,'on')
                obj.Back.EdgeColor = obj.PrvtBackEdgeColor;
            else
                obj.Back.EdgeColor = obj.OffBackEdgeColor;            
            end
            
        end
        function Value = get.BackEdgeColor(obj)
            Value = obj.PrvtBackEdgeColor;
        end
        
        function set.BackLineWidth(obj,Value)

            % Checking input values and updating private property
            if Value(1) <= 0
                error('Line width value must be greater than 0');
            end
            obj.PrvtBackLineWidth = Value(1);

            % Updating relevant gage components
            obj.Back.LineWidth = obj.PrvtBackLineWidth;
            
        end
        function Value = get.BackLineWidth(obj)
            Value = obj.PrvtBackLineWidth;
        end

        function set.BackVisible(obj,Value)
            
            % Checking input values and updating private property
            if ~strcmpi(Value,'on') &&  ~strcmpi(Value,'off')
                error('''BackVisible'' property must be ''on'' or ''off''');
            end
            obj.PrvtBackVisible = Value;
            
            % Updating relevant gage components
            obj.Back.Visible = obj.PrvtBackVisible;
            
        end
        function Value = get.BackVisible(obj)
            Value = obj.PrvtBackVisible;
        end
        
        
% SET/GET METHODS FOR GAGE COLOR BAND

        function set.BandValues(obj,Value)
            
            % Checking input values
            if numel(Value) == 1
                error('''BandValues'' must have at least two elements' );
            end
            if any(diff(Value)<=0)
                error('''BandValues'' must be in ascending order');
            end
            if size(Value,1) > size(Value,2)
                Value = Value';
            end
            
            % Checking for change in number of band values
            if length(Value) ~= length(obj.PrvtBandValues)
                Option = 'new';
            else
                Option = 'updated';
            end
            
            % Updating private property
            obj.PrvtBandValues = Value;
            
            % Updating BandColor property
            if size(obj.BandColor,1) > length(obj.PrvtBandValues)
                obj.PrvtBandColor = obj.BandColor(1:length(obj.PrvtBandValues),:);
            elseif size(obj.BandColor,1) < length(obj.PrvtBandValues)
                obj.PrvtBandColor = [obj.BandColor; ones(length(obj.PrvtBandValues)-size(obj.BandColor,1),3)];
            end
            
            % Updating 'off' band color
            obj.OffBandColor = 0.85*ones(size(obj.PrvtBandColor));
            
            % Calculating angle arrays
           [obj.ThtDial,obj.ThtBand,obj.ThtTick,obj.ThtMinorTick] = obj.calculateAngles(...
               obj.DialAngle,obj.BandValues,obj.DialValues,obj.DialMinorValues);
           
            % Updating relevant gage components
            obj.drawBand(Option)
            obj.drawDial(Option)
            obj.drawMinorTick(Option)
            obj.drawTick(Option)
            obj.drawLabel(Option)
            obj.drawNeedle
            
        end
        function Value = get.BandValues(obj)
            Value = obj.PrvtBandValues;
        end
        
        function set.BandColor(obj,Value)
            
            % Checking input values and updating private property
            if size(Value,1) ~= length(obj.BandValues)  ||  size(Value,2) ~= 3
                error('Size of ''BandValues'' array must be (length(BandValues),3)');
            end
            obj.PrvtBandColor = Value;
            
            % Updating relevant gage components
            obj.drawBand('updated')
            
        end   
        function Value = get.BandColor(obj)
            Value = obj.PrvtBandColor;
        end
        
        function set.BandPosition(obj,Value)
            
            if numel(Value) ~= 2
                error('''BandPosition'' property must have 2 values');
            end
            if size(Value,1) > size(Value,2)
                Value = Value';
            end
            obj.PrvtBandPosition = Value;
            
            % Updating relevant gage components
            obj.drawBand('updated')
            
        end   
        function Value = get.BandPosition(obj)
            Value = obj.PrvtBandPosition;
        end
        
        function set.BandVisible(obj,Value)
            
            % Checking input values and updating private property
            if ~strcmpi(Value,'on') &&  ~strcmpi(Value,'off')
                error('''BandVisible'' property must be ''on'' or ''off''');
            end
            obj.PrvtBandVisible = Value;
            
            % Updating relevant gage components
            for i = 1:length(obj.Band)
                obj.Band{i}.Visible = obj.PrvtBandVisible;
            end
            
        end
        function Value = get.BandVisible(obj)
            Value = obj.PrvtBandVisible;
        end
        
        
% SET/GET METHODS FOR GAGE DIAL

        function set.DialAngle(obj,Value)
            
            % Checking input values and updating private property
            if numel(Value) ~= 2
                error('''DialAngle'' property must have 2 values');
            end
            if size(Value,1) > size(Value,2)
                Value = Value';
            end
            if Value(1) >= Value(2)
                error('Second angle value must be greater than first angle value');
            end
            obj.PrvtDialAngle = Value;
            
            % Calculating angle arrays
           [obj.ThtDial,obj.ThtBand,obj.ThtTick,obj.ThtMinorTick] = obj.calculateAngles(...
               obj.DialAngle,obj.BandValues,obj.DialValues,obj.DialMinorValues);
           
            % Updating needle position
            obj.NeedlePosition = obj.DialValues(1);
            
            % Updating relevant gage components
            obj.drawBand('new')
            obj.drawDial('new')
            obj.drawMinorTick('new')
            obj.drawTick('new')
            obj.drawLabel('new')
            obj.drawNeedle
            
        end
        function Value = get.DialAngle(obj)
            Value = obj.PrvtDialAngle;
        end
        
        function set.DialValues(obj,Value)
            
            % Checking input values
            if numel(Value) == 1
                error('''DialValues'' must have at least two elements' );
            end
            if any(diff(Value)<=0)
                error('''DialValues'' must be in ascending order');
            end
            if size(Value,1) > size(Value,2)
                Value = Value';
            end
            
            % Checking for change in number of dial values
            if length(Value) ~= length(obj.PrvtDialValues)
                Option = 'new';
            else
                Option = 'updated';
            end
            
            % Updating private property
            obj.PrvtDialValues = Value;
            
            % Calculating new color band valules
            BandValueNew = zeros(1,length(obj.BandValues));
            BandValueNew(1) = obj.PrvtDialValues(1) + (obj.PrvtDialValues(end)-obj.PrvtDialValues(1))/...
                (obj.ThtTick(end)-obj.ThtTick(1)) * (obj.ThtBand{1}(1)-obj.ThtTick(1));
            for i = 1:length(obj.ThtBand)
                BandValueNew(i+1) = obj.PrvtDialValues(1) + (obj.PrvtDialValues(end)-obj.PrvtDialValues(1))/...
                    (obj.ThtTick(end)-obj.ThtTick(1)) * (obj.ThtBand{i}(end)-obj.ThtTick(1));
            end
            
            % Treating funky case where band has repeat values
            % due to very small separation between band values
            k = find(diff(BandValueNew)==0);
            for i = 1:length(k)
                BandValueNew(k(i)) = BandValueNew(k(i)+1) - (BandValueNew(end)-BandValueNew(1))/100;
            end
            
            % Updating minor tick values
            obj.DialMinorValues = obj.PrvtDialValues(1): ...
                (obj.PrvtDialValues(2)-obj.PrvtDialValues(1))/5: ...
                 obj.PrvtDialValues(end);
            
             % Calculating angle arrays
             [obj.ThtDial,obj.ThtBand,obj.ThtTick,obj.ThtMinorTick] = obj.calculateAngles(...
                 obj.DialAngle,BandValueNew,obj.DialValues,obj.DialMinorValues);
             
            % Updating needle position
            obj.NeedlePosition = obj.DialValues(1);
            
            % Updating relevant gage components
            obj.drawBand(Option)
            obj.drawDial(Option)
            obj.drawMinorTick(Option)
            obj.drawTick(Option)
            obj.drawLabel('new')
            obj.drawDisplay('updated')
            obj.drawNeedle
            
            % Updating color band values
            obj.BandValues = BandValueNew;
            
        end
        function Value = get.DialValues(obj)
            Value = obj.PrvtDialValues;
        end
        
        function set.DialMinorValues(obj,Value)
            
            % Checking input values and updating private property
            if numel(Value) == 1
                error('''DialMinorValues'' must have at least two elements' );
            end
            if any(diff(Value)<=0)
                error('''DialMinorValues'' must be in ascending order');
            end
            if size(Value,1) > size(Value,2)
                Value = Value';
            end
            
            % Checking for change in number of dial values
            if length(Value) ~= length(obj.PrvtDialMinorValues)
                Option = 'new';
            else
                Option = 'updated';
            end
            
            % Updating private property
            obj.PrvtDialMinorValues = Value;
            
            % Calculating angle arrays
           [obj.ThtDial,obj.ThtBand,obj.ThtTick,obj.ThtMinorTick] = obj.calculateAngles(...
               obj.DialAngle,obj.BandValues,obj.DialValues,obj.DialMinorValues);
           
            % Updating relevant gage components
            obj.drawMinorTick(Option)
            
        end
        function Value = get.DialMinorValues(obj)
            Value = obj.PrvtDialMinorValues;
        end
        
        function set.DialColor(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtDialColor = Value;
            
            % Updating relevant gage components
            if strcmpi(obj.Enable,'on')
                obj.Dial.Color = obj.PrvtDialColor;
            else
                obj.Dial.Color = obj.OffDialColor;
            end
            
        end
        function Value = get.DialColor(obj)
            Value = obj.PrvtDialColor;
        end
        
        function set.DialLineWidth(obj,Value)

            % Checking input values and updating private property
            if Value(1) <= 0
                error('''DialLineWidth'' width value must be greater than 0');
            end
            obj.PrvtDialLineWidth = Value(1);
            
            % Updating relevant gage components
            obj.Dial.LineWidth = obj.PrvtDialLineWidth;
            
        end
        function Value = get.DialLineWidth(obj)
            Value = obj.PrvtDialLineWidth;
        end
        
        function set.DialVisible(obj,Value)
            
            % Checking input values and updating private property
            if ~strcmpi(Value,'on') &&  ~strcmpi(Value,'off')
                error('''DialVisible'' property must be ''on'' or ''off''');
            end
            obj.PrvtDialVisible = Value;
            
            % Updating relevant gage components
            obj.Dial.Visible = obj.PrvtDialVisible;
            
        end
        function Value = get.DialVisible(obj)
            Value = obj.PrvtDialVisible;
        end
        
        
% SET/GET METHODS FOR GAGE DIAL LABELS

        function set.LabelString(obj,Value)
            
            % Checking input values and updating private property
            if length(Value) ~= length(obj.DialValues)
                error('''LabelString'' property must have the same length as ''DialValues''');
            end
            obj.PrvtLabelString = Value;
            
            % Updating relevant gage components
            for i = 1:length(obj.Label)
                obj.Label{i}.String = obj.PrvtLabelString{i};
            end
            
        end
        function Value = get.LabelString(obj)
            Value = obj.PrvtLabelString;
        end
        
        function set.LabelFontName(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtLabelFontName = Value;
            
            % Updating relevant gage components
            for i = 1:length(obj.Label)
                obj.Label{i}.FontName = obj.PrvtLabelFontName;
            end
            
        end
        function Value = get.LabelFontName(obj)
            Value = obj.PrvtLabelFontName;
        end
        
        function set.LabelFontSize(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtLabelFontSize = Value;
            
            % Updating relevant gage components
            for i = 1:length(obj.Label)
                obj.Label{i}.FontSize = obj.PrvtLabelFontSize;
            end
            
            % Updating relevant gage components
            obj.drawLabel('updated');
            obj.drawNeedle
            
        end
        function Value = get.LabelFontSize(obj)
            Value = obj.PrvtLabelFontSize;
        end
        
        function set.LabelFontWeight(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtLabelFontWeight = Value;
            
            % Updating relevant gage components
            for i = 1:length(obj.Label)
                obj.Label{i}.FontWeight = obj.PrvtLabelFontWeight;
            end
            
        end
        function Value = get.LabelFontWeight(obj)
            Value = obj.PrvtLabelFontWeight;
        end
        
        function set.LabelColor(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtLabelColor = Value;
            
            % Updating relevant gage components
            if strcmpi(obj.Enable,'on')
                for i = 1:length(obj.Label)
                    obj.Label{i}.Color = obj.PrvtLabelColor;
                end
            else
                for i = 1:length(obj.Label)
                    obj.Label{i}.Color = obj.OffLabelColor;
                end
            end
            
        end
        function Value = get.LabelColor(obj)
            Value = obj.PrvtLabelColor;
        end
        
        function set.LabelDir(obj,Value)
            
            % Checking input values and updating private property
            if ~strcmpi(Value,'in') &&  ~strcmpi(Value,'out')  &&  ...
                    ~strcmpi(Value,'center')
                error('''LabelDir'' property must be ''in'' , ''out'' or ''center''');
            end
            obj.PrvtLabelDir = Value;
            
            % Updating relevant gage components
            obj.drawLabel('updated');
            obj.drawNeedle
            
        end
        function Value = get.LabelDir(obj)
            Value = obj.PrvtLabelDir;
        end
        
        function set.LabelVisible(obj,Value)
            
            % Checking input values and updating private property
            if ~strcmpi(Value,'on') &&  ~strcmpi(Value,'off')
                error('''LabelVisible'' property must be ''on'' or ''off''');
            end
            obj.PrvtLabelVisible = Value;
            
            % Updating relevant gage components
            for i = 1:length(obj.Label)
                obj.Label{i}.Visible = obj.PrvtLabelVisible;
            end
            
        end
        function Value = get.LabelVisible(obj)
            Value = obj.PrvtLabelVisible;
        end
        
        function set.LabelOrientation(obj,Value)
            
            % Checking input values and updating private property
            if ~strcmpi(Value,'horizontal') &&  ~strcmpi(Value,'radial')  &&  ...
                    ~strcmpi(Value,'tangential')
                error('''LabelOrientation'' property must be ''horizontal'' , ''radial'' or ''tangential''');
            end
            obj.PrvtLabelOrientation = Value;
            
            % Calculating label orientation
            LabelRotation = zeros(1,length(obj.ThtTick));
            switch obj.PrvtLabelOrientation
                case 'horizontal'

                case 'radial'
                    for i = 1:length(obj.ThtTick)
                        LabelRotation(i) = 180/pi*obj.ThtTick(i)-90;
                    end
                case 'tangential'
                    for i = 1:length(obj.ThtTick)
                        LabelRotation(i) = 180/pi*obj.ThtTick(i);
                    end
            end
            
            % Updating relevant gage components
            for i = 1:length(obj.Label)
                obj.Label{i}.Rotation = LabelRotation(i);
            end
            
        end
        function Value = get.LabelOrientation(obj)
            Value = obj.PrvtLabelOrientation;
        end
        
        
% SET/GET METHODS FOR GAGE MAJOR TICKS

        function set.TickColor(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtTickColor = Value;
            
            % Updating relevant gage components
            if strcmpi(obj.Enable,'on')
                for i = 1:length(obj.Tick)
                    obj.Tick{i}.Color = obj.PrvtTickColor;
                end
            else
                for i = 1:length(obj.Tick)
                    obj.Tick{i}.Color = obj.OffTickColor;
                end
            end
            
        end
        function Value = get.TickColor(obj)
            Value = obj.PrvtTickColor;
        end
        
        function set.TickDir(obj,Value)
            
            % Checking input values and updating private property
            if ~strcmpi(Value,'in') &&  ~strcmpi(Value,'out')  &&  ...
                    ~strcmpi(Value,'center')
                error('''TickDir'' property must be ''in'' , ''out'' or ''center''');
            end
            obj.PrvtTickDir = Value;
            
            % Updating relevant gage components
            obj.drawTick('updated')
            obj.drawLabel('updated')
            
        end
        function Value = get.TickDir(obj)
            Value = obj.PrvtTickDir;
        end
        
        function set.TickLength(obj,Value)

            % Checking input values and updating private property
            if Value(1) <= 0  ||  Value(1) > 1
                error('''TickLength'' property must be between 0 and 1');
            end
            obj.PrvtTickLength = Value(1);
            
            % Updating relevant gage components
            obj.drawTick('updated')
            obj.drawLabel('updated')
            
        end
        function Value = get.TickLength(obj)
            Value = obj.PrvtTickLength;
        end
        
        function set.TickLineWidth(obj,Value)

            % Checking input values and updating private property
            if Value(1) <= 0
                error('''TickLineWidth'' value must be greater than 0');
            end
            obj.PrvtTickLineWidth = Value(1);
            
            % Updating relevant gage components
            for i = 1:length(obj.Tick)
                obj.Tick{i}.LineWidth = obj.PrvtTickLineWidth;
            end
            
        end
        function Value = get.TickLineWidth(obj)
            Value = obj.PrvtTickLineWidth;
        end
        
        function set.TickPosition(obj,Value)

            % Checking input values and updating private property
            if Value(1) <= 0
                error('''TickPosition'' property must be greater than 0');
            end
            obj.PrvtTickPosition = Value(1);
            
            % Updating relevant gage components
            obj.drawTick('updated')
            obj.drawLabel('updated')
            
        end
        function Value = get.TickPosition(obj)
            Value = obj.PrvtTickPosition;
        end
        
        function set.TickVisible(obj,Value)
            
            % Checking input values and updating private property
            if ~strcmpi(Value,'on') &&  ~strcmpi(Value,'off')
                error('''TickVisible'' property must be ''on'' or ''off''');
            end
            obj.PrvtTickVisible = Value;
            
            % Updating relevant gage components
            for i = 1:length(obj.Tick)
                obj.Tick{i}.Visible = obj.PrvtTickVisible;
            end
            
        end
        function Value = get.TickVisible(obj)
            Value = obj.PrvtTickVisible;
        end
        
        
% SET/GET METHODS FOR GAGE MINOR TICKS

        function set.MinorTickColor(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtMinorTickColor = Value;
            
            % Updating relevant gage components
            if strcmpi(obj.Enable,'on')
                for i = 1:length(obj.MinorTick)
                    obj.MinorTick{i}.Color = obj.PrvtMinorTickColor;
                end
            else
                for i = 1:length(obj.MinorTick)
                    obj.MinorTick{i}.Color = obj.OffMinorTickColor;
                end
            end
            
        end
        function Value = get.MinorTickColor(obj)
            Value = obj.PrvtMinorTickColor;
        end
        
        function set.MinorTickDir(obj,Value)
            
            % Checking input values and updating private property
            if ~strcmpi(Value,'in') &&  ~strcmpi(Value,'out')  &&  ...
                    ~strcmpi(Value,'center')
                error('''MinorTickDir'' property must be ''in'' , ''out'' or ''center''');
            end
            obj.PrvtMinorTickDir = Value;
            
            % Updating relevant gage components
            obj.drawMinorTick('updated')
            
        end
        function Value = get.MinorTickDir(obj)
            Value = obj.PrvtMinorTickDir;
        end
        
        function set.MinorTickLength(obj,Value)

            % Checking input values and updating private property
            if Value(1) <= 0  ||  Value(1) > 1
                error('''MinorTickLength'' property must be between 0 and 1');
            end
            obj.PrvtMinorTickLength = Value(1);
            
            % Updating relevant gage components
            obj.drawMinorTick('updated')
            
        end
        function Value = get.MinorTickLength(obj)
            Value = obj.PrvtMinorTickLength;
        end
        
        function set.MinorTickLineWidth(obj,Value)

            % Checking input values and updating private property
            if Value(1) <= 0
                error('''MinorTickLineWidth'' width value must be greater than 0');
            end
            obj.PrvtMinorTickLineWidth = Value(1);
            
            % Updating relevant gage components
            for i = 1:length(obj.MinorTick)
                obj.MinorTick{i}.LineWidth = obj.PrvtMinorTickLineWidth;
            end
            
        end
        function Value = get.MinorTickLineWidth(obj)
            Value = obj.PrvtMinorTickLineWidth;
        end
        
        function set.MinorTickPosition(obj,Value)

            % Checking input values and updating private property
            if Value(1) <= 0
                error('''MinorTickPosition'' property must be greater than 0');
            end
            obj.PrvtMinorTickPosition = Value(1);
            
            % Updating relevant gage components
            obj.drawMinorTick('updated')
            
        end
        function Value = get.MinorTickPosition(obj)
            Value = obj.PrvtMinorTickPosition;
        end
        
        function set.MinorTickVisible(obj,Value)
            
            % Checking input values and updating private property
            if ~strcmpi(Value,'on') &&  ~strcmpi(Value,'off')
                error('''MinorTickVisible'' property must be ''on'' or ''off''');
            end
            obj.PrvtMinorTickVisible = Value;
            
            % Updating relevant gage components
            for i = 1:length(obj.MinorTick)
                obj.MinorTick{i}.Visible = obj.PrvtMinorTickVisible;
            end
            
        end
        function Value = get.MinorTickVisible(obj)
            Value = obj.PrvtMinorTickVisible;
        end
        
        
% SET/GET METHODS FOR GAGE KNOB

        function set.KnobSize(obj,Value)
            
            % Checking input values and updating private property
            if Value(1) <= 0
                error('''KnobSize'' property must be greater than 0');
            end
            obj.PrvtKnobSize = Value(1);
            
            % Updating relevant gage components
            obj.drawKnob('updated')
            
        end
        function Value = get.KnobSize(obj)
            Value = obj.PrvtKnobSize;
        end
        
        function set.KnobColor(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtKnobColor = Value;
            
            % Updating relevant gage components
            if strcmpi(obj.Enable,'on')
                obj.Knob.FaceColor = obj.PrvtKnobColor;
            else
                obj.Knob.FaceColor = obj.OffKnobColor;
            end
            
        end
        function Value = get.KnobColor(obj)
            Value = obj.PrvtKnobColor;
        end
        
        function set.KnobEdgeColor(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtKnobEdgeColor = Value;
            
            % Updating relevant gage components
            if strcmpi(obj.Enable,'on')
                obj.Knob.EdgeColor = obj.PrvtKnobEdgeColor;
            else
                obj.Knob.EdgeColor = obj.OffKnobEdgeColor;  
            end
            
        end
        function Value = get.KnobEdgeColor(obj)
            Value = obj.PrvtKnobEdgeColor;
        end
        
        function set.KnobLineWidth(obj,Value)

            % Checking input values and updating private property
            if Value(1) <= 0
                error('''KnobLineWidth'' width value must be greater than 0');
            end
            obj.PrvtKnobLineWidth = Value(1);
            
            % Updating relevant gage components
            obj.Knob.LineWidth = obj.PrvtKnobLineWidth;
            
        end
        function Value = get.KnobLineWidth(obj)
            Value = obj.PrvtKnobLineWidth;
        end
        
        
% SET/GET METHODS FOR GAGE NEEDLE

        function set.NeedleColor(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtNeedleColor = Value;
            
            % Updating relevant gage components
            obj.Needle.FaceColor = obj.PrvtNeedleColor;
            
        end
        function Value = get.NeedleColor(obj)
            Value = obj.PrvtNeedleColor;
        end
        
        function set.NeedleEdgeColor(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtNeedleEdgeColor = Value;
            
            % Updating relevant gage components
            obj.Needle.EdgeColor = obj.PrvtNeedleEdgeColor;
            
        end
        function Value = get.NeedleEdgeColor(obj)
            Value = obj.PrvtNeedleEdgeColor;
        end
        
        function set.NeedleLineWidth(obj,Value)

            % Checking input values and updating private property
            if Value(1) <= 0
                error('''NeedleLineWidth'' width value must be greater than 0');
            end
            obj.PrvtNeedleLineWidth = Value(1);
            
            % Updating relevant gage components
            obj.Needle.LineWidth = obj.PrvtNeedleLineWidth;
            
        end
        function Value = get.NeedleLineWidth(obj)
            Value = obj.PrvtNeedleLineWidth;
        end
        
        function set.NeedleLength(obj,Value)

            % Checking input values and updating private property
            if Value(1) <= 0
                error('''NeedleLength'' property must be greater than 0');
            end
            obj.PrvtNeedleLength = Value(1);
            
            % Clearing needle points property
            obj.NeedlePoints = [];
            
            % Updating relevant gage components
            obj.drawNeedle;
            
        end
        function Value = get.NeedleLength(obj)
            Value = obj.PrvtNeedleLength;
        end
        
        function set.NeedleWidth(obj,Value)

            % Checking input values and updating private property
            if numel(Value) < 1  ||  numel(Value) > 2
                error('''NeedleWidth'' property must have 1 or 2 values');
            end
            if any(Value<=0)  ||  any(Value>1)
                error('''NeedleWidth'' property must be between 0 and 1');
            end
            if size(Value,1) > size(Value,2)
                Value = Value';
            end
            obj.PrvtNeedleWidth = Value;
            
            % Clearing needle points property
            obj.NeedlePoints = [];
            
            % Updating relevant gage components
            obj.drawNeedle;
            
        end
        function Value = get.NeedleWidth(obj)
            Value = obj.PrvtNeedleWidth;
        end
        
        function set.NeedleOffset(obj,Value)

            % Checking input values and updating private property
            if Value(1) <= -1  ||  Value(1) > 1
                error('''NeedleOffset'' property must be between -1 and 1');
            end
            obj.PrvtNeedleOffset = Value(1);
            
            % Updating relevant gage components
            obj.drawNeedle;
            
        end
        function Value = get.NeedleOffset(obj)
            Value = obj.PrvtNeedleOffset;
        end

        function set.NeedlePoints(obj,Value)

            % Checking input values and updating private property
            if ~isempty(Value)
                if size(Value,1) < 3  ||  size(Value,2) ~= 2
                    error('''NeedlePoints'' property must have at least 3 points and be size (N,2)');
                end
            end
            obj.PrvtNeedlePoints = Value;
            
            % Updating relevant gage components
            obj.drawNeedle;
            
        end
        function Value = get.NeedlePoints(obj)
            Value = obj.PrvtNeedlePoints;
        end
        
        function set.NeedlePosition(obj,Value)
           
            % Checking value limits
            if Value < min(obj.DialMinorValues)
                obj.PrvtNeedlePosition = min(obj.DialMinorValues);
            elseif Value > max(obj.DialMinorValues)
                obj.PrvtNeedlePosition = max(obj.DialMinorValues);
            else
                obj.PrvtNeedlePosition = Value;
            end
            
            % Updating numerical display
            obj.drawDisplay('value')
            
            % Moving needle
            obj.moveNeedle(obj.PrvtNeedlePosition)
            
        end
        function Value = get.NeedlePosition(obj)
            Value = obj.PrvtNeedlePosition;
        end
        
        
% SET/GET METHODS FOR NUMERICAL DISPLAY

        function set.DisplayFontName(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtDisplayFontName = Value;
            
            % Updating relevant gage components
            obj.Display{2}.FontName = obj.PrvtDisplayFontName;
            
        end
        function Value = get.DisplayFontName(obj)
            Value = obj.PrvtDisplayFontName;
        end
        
        function set.DisplayFontSize(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtDisplayFontSize = Value;
            
            % Updating relevant gage components
            obj.Display{2}.FontSize = obj.PrvtDisplayFontSize;
            
        end
        function Value = get.DisplayFontSize(obj)
            Value = obj.PrvtDisplayFontSize;
        end
        
        function set.DisplayFontWeight(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtDisplayFontWeight = Value;
            
            % Updating relevant gage components
            obj.Display{2}.FontWeight = obj.PrvtDisplayFontWeight;
            
        end
        function Value = get.DisplayFontWeight(obj)
            Value = obj.PrvtDisplayFontWeight;
        end
        
        function set.DisplayForegroundColor(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtDisplayForegroundColor = Value;
            
            % Updating relevant gage components
            if strcmpi(obj.Enable,'on')
                obj.Display{2}.Color = obj.PrvtDisplayForegroundColor;
            else
                obj.Display{2}.Color = obj.OffDisplayForegroundColor;
            end
            
        end
        function Value = get.DisplayForegroundColor(obj)
            Value = obj.PrvtDisplayForegroundColor;
        end
        
        function set.DisplayBackgroundColor(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtDisplayBackgroundColor = Value;
            
            % Updating relevant gage components
            if strcmpi(obj.Enable,'on')
                obj.Display{1}.FaceColor = obj.PrvtDisplayBackgroundColor;
            else
                obj.Display{1}.FaceColor = obj.OffDisplayBackgroundColor;
            end
        end
        function Value = get.DisplayBackgroundColor(obj)
            Value = obj.PrvtDisplayBackgroundColor;
        end
        
        function set.DisplayEdgeColor(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtDisplayEdgeColor = Value;
            
            % Updating relevant gage components
            if strcmpi(obj.Enable,'on')
                obj.Display{1}.EdgeColor = obj.PrvtDisplayEdgeColor;
            else
                obj.Display{1}.EdgeColor = obj.OffDisplayEdgeColor;
            end
        end
        function Value = get.DisplayEdgeColor(obj)
            Value = obj.PrvtDisplayEdgeColor;
        end
        
        function set.DisplayLineWidth(obj,Value)

            % Checking input values and updating private property
            if Value(1) <= 0
                error('''DisplayLineWidth'' width value must be greater than 0');
            end
            obj.PrvtDisplayLineWidth = Value(1);
            
            % Updating relevant gage components
            obj.Display{1}.LineWidth = obj.PrvtDisplayLineWidth;
            
        end
        function Value = get.DisplayLineWidth(obj)
            Value = obj.PrvtDisplayLineWidth;
        end
        
        function set.DisplayPosition(obj,Value)
           
            % Checking value limits
            % Checking input values and updating private property
            if size(Value,1) ~= 1  ||  size(Value,2) ~= 4
                error('Size of position array must be (1,4)');
            end
            obj.PrvtDisplayPosition = Value;
            
            % Updating relevant gage components
            obj.drawDisplay('updated')
            
        end
        function Value = get.DisplayPosition(obj)
            Value = obj.PrvtDisplayPosition;
        end
        
        function set.DisplayVisible(obj,Value)
            
            % Checking input values and updating private property
            if ~strcmpi(Value,'on') &&  ~strcmpi(Value,'off')
                error('''DisplayVisible'' property must be ''on'' or ''off''');
            end
            obj.PrvtDisplayVisible = Value;
            
            % Updating relevant gage components
            obj.Display{1}.Visible = obj.PrvtDisplayVisible;
            obj.Display{2}.Visible = obj.PrvtDisplayVisible;
            
        end
        function Value = get.DisplayVisible(obj)
            Value = obj.PrvtDisplayVisible;
        end
        
        
% SET/GET METHODS FOR UNIT DISPLAY

        function set.UnitString(obj,Value)
           
            % Checking value limits
            obj.PrvtUnitString = Value;
            
            % Updating relevant gage components
            obj.Unit.String = obj.PrvtUnitString;
            
        end
        function Value = get.UnitString(obj)
            Value = obj.PrvtUnitString;
        end
        
        function set.UnitFontName(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtUnitFontName = Value;
            
            % Updating relevant gage components
            obj.Unit.FontName = obj.PrvtUnitFontName;
            
        end
        function Value = get.UnitFontName(obj)
            Value = obj.PrvtUnitFontName;
        end
        
        function set.UnitFontSize(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtUnitFontSize = Value;
            
            % Updating relevant gage components
            obj.Unit.FontSize = obj.PrvtUnitFontSize;
            
        end
        function Value = get.UnitFontSize(obj)
            Value = obj.PrvtUnitFontSize;
        end
        
        function set.UnitFontWeight(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtUnitFontWeight = Value;
            
            % Updating relevant gage components
            obj.Unit.FontWeight = obj.PrvtUnitFontWeight;
            
        end
        function Value = get.UnitFontWeight(obj)
            Value = obj.PrvtUnitFontWeight;
        end
        
        function set.UnitColor(obj,Value)
            
            % Checking input values and updating private property
            obj.PrvtUnitColor = Value;
            
            % Updating relevant gage components
            if strcmpi(obj.Enable,'on')
                obj.Unit.Color = obj.PrvtUnitColor;
            else
                obj.Unit.Color = obj.OffUnitColor;
            end
            
        end
        function Value = get.UnitColor(obj)
            Value = obj.PrvtUnitColor;
        end
        
        function set.UnitPosition(obj,Value)
           
            % Checking value limits
            if size(Value,1) ~= 1  ||  size(Value,2) ~= 2
                error('Size of position array must be (1,2)');
            end
            obj.PrvtUnitPosition = Value;
            
            % Updating relevant gage components
            obj.drawUnit('updated')
            
        end
        function Value = get.UnitPosition(obj)
            Value = obj.PrvtUnitPosition;
        end
        
        function set.UnitVisible(obj,Value)
            
            % Checking input values and updating private property
            if ~strcmpi(Value,'on') &&  ~strcmpi(Value,'off')
                error('''UnitVisible'' property must be ''on'' or ''off''');
            end
            obj.PrvtUnitVisible = Value;
            
            % Updating relevant gage components
            obj.Unit.Visible = obj.PrvtUnitVisible;
            
        end
        function Value = get.UnitVisible(obj)
            Value = obj.PrvtUnitVisible;
        end

    end
    
    
% =========================================================================           
%                  DEFINING GAGE COMPONENTS UPDATE METHODS    
% =========================================================================           
    
    methods (Hidden = true)
        
        function drawBack(obj,varargin)
           
            % Checking for update option
            Option = 'new';
            if nargin == 2
                if strcmpi(varargin{1},'new')
                    Option = 'new';
                elseif strcmpi(varargin{1},'updated')
                    Option = 'updated';
                end
            end

            % Creating auxiliary angle
            ThtAux = 0:pi/180:2*pi;
            
            % Calculating back panel points
            PtBack = [cos(ThtAux); sin(ThtAux)];

            % Creating auxiliary color
            if strcmpi(obj.Enable,'on')
                ColorAux1 = obj.BackColor;
                ColorAux2 = obj.BackEdgeColor;
            else
                ColorAux1 = obj.OffBackColor;
                ColorAux2 = obj.OffBackEdgeColor;
            end

            % Checking for existing back panel
            if strcmp(Option,'new')
                axes(obj.Axis)
                obj.Back = patch(PtBack(1,:),PtBack(2,:),ColorAux1);
                set(obj.Back,...
                    'FaceLighting','flat',...
                    'EdgeColor',ColorAux2,...
                    'LineWidth',obj.BackLineWidth);
            else
                obj.Back.XData = PtBack(1,:);
                obj.Back.YData = PtBack(2,:);
            end
            
        end
        
        function drawBand(obj,varargin)
            
            % Checking for update option
            Option = 'new';
            if nargin == 2
                if strcmpi(varargin{1},'new')
                    Option = 'new';
                elseif strcmpi(varargin{1},'updated')
                    Option = 'updated';
                end
            end

            % Creating auxiliary color
            if strcmpi(obj.Enable,'on')
                ColorAux = obj.BandColor;
            else
                ColorAux = obj.OffBandColor;
            end

            % Looping through color band segments
            PtBand = cell(1,length(obj.ThtBand));
            PtColor = cell(1,length(obj.ThtBand));
            for i = 1:length(PtBand)

                % Calculating band points
                PtBandBot = obj.BandPosition(1)*[cos(obj.ThtBand{i}); sin(obj.ThtBand{i})];
                PtBandTop = obj.BandPosition(2)*[cos(obj.ThtBand{i}); sin(obj.ThtBand{i})];
                PtBand{i} = [PtBandBot fliplr(PtBandTop)];
                
                % Calculating band vertex colors
                PtColorBot = ones(length(obj.ThtBand{i}),1)*ColorAux(i,:) + ...
                    (0:length(obj.ThtBand{i})-1)'* ...
                    (ColorAux(i+1,:)-ColorAux(i,:))/(length(obj.ThtBand{i})-1);
                PtColorTop = flipud(PtColorBot);
                PtColor{i} = [PtColorBot; PtColorTop];
                
            end
            
            % Checking for change in number of ticks
            if strcmp(Option,'new')
                
                % Deleting old ticks
                for i = 1:length(obj.Band)
                    delete(obj.Band{i})
                end
                obj.Band = cell(1,length(obj.ThtBand));
                
                % Creating new ticks
                axes(obj.Axis)
                for i = 1:length(obj.Band)
                    obj.Band{i} = patch('XData',PtBand{i}(1,:),'YData',PtBand{i}(2,:));
                    set(obj.Band{i},...
                        'FaceVertexCData',PtColor{i},...
                        'FaceColor','interp',...
                        'FaceLighting','flat',...
                        'EdgeColor','none');
                end
                
            else
                
                % Updating existing ticks
                for i = 1:length(obj.Band)
                    obj.Band{i}.XData = PtBand{i}(1,:);
                    obj.Band{i}.YData = PtBand{i}(2,:);
                    obj.Band{i}.FaceVertexCData = PtColor{i};
                end
                
            end
            
        end
        
        function drawKnob(obj,varargin)
            
            % Checking for update option
            Option = 'new';
            if nargin == 2
                if strcmpi(varargin{1},'new')
                    Option = 'new';
                elseif strcmpi(varargin{1},'updated')
                    Option = 'updated';
                end
            end

            % Creating auxiliary angle
            ThtAux = 0:pi/180:2*pi;
            
            % Calculating knob points
            PtKnob = 0.5*obj.KnobSize*[cos(ThtAux); sin(ThtAux)];
            
            % Creating auxiliary color
            if strcmpi(obj.Enable,'on')
                ColorAux1 = obj.KnobColor;
                ColorAux2 = obj.KnobEdgeColor;
            else
                ColorAux1 = obj.OffKnobColor;
                ColorAux2 = obj.OffKnobEdgeColor;
            end

            % Checking for existing knob
            if strcmpi(Option,'new')
                
                axes(obj.Axis)
                obj.Knob = patch(PtKnob(1,:),PtKnob(2,:),ColorAux1);
                set(obj.Knob,...
                    'FaceLighting','flat',...
                    'EdgeColor',ColorAux2,...
                    'LineWidth',obj.KnobLineWidth);
                
                % Updating needle
                obj.drawNeedle

            else
                obj.Knob.XData = PtKnob(1,:);
                obj.Knob.YData = PtKnob(2,:);
            end
            
        end            
        
        function drawDial(obj,varargin)
            
            % Checking for update option
            Option = 'new';
            if nargin == 2
                if strcmpi(varargin{1},'new')
                    Option = 'new';
                elseif strcmpi(varargin{1},'updated')
                    Option = 'updated';
                end
            end
            
            % Calculating dial points angles
            PtDial = [cos(obj.ThtDial); sin(obj.ThtDial)];
            
            % Creating auxiliary color
            if strcmpi(obj.Enable,'on')
                ColorAux = obj.DialColor;
            else
                ColorAux = obj.OffDialColor;
            end

            % Checking for existing dial
            if strcmpi(Option,'new')
                
                if ~isempty(obj.Dial)
                    delete(obj.Dial)
                end
                
                obj.Dial = plot(obj.Axis,PtDial(1,:),PtDial(2,:));
                set(obj.Dial,...
                    'Color',ColorAux,...
                    'LineWidth',obj.DialLineWidth);
            else
                obj.Dial.XData = PtDial(1,:);
                obj.Dial.YData = PtDial(2,:);
            end

        end
        
        function drawMinorTick(obj,varargin)
            
            % Checking for update option
            Option = 'new';
            if nargin == 2
                if strcmpi(varargin{1},'new')
                    Option = 'new';
                elseif strcmpi(varargin{1},'updated')
                    Option = 'updated';
                end
            end
            
            % Creating tick end point radius
            if strcmpi(obj.MinorTickDir,'in')
                R1 = obj.MinorTickPosition-obj.MinorTickLength;
                R2 = obj.MinorTickPosition;
            elseif strcmpi(obj.MinorTickDir,'out')
                R1 = obj.MinorTickPosition;
                R2 = obj.MinorTickPosition+obj.MinorTickLength;
            else
                R1 = obj.MinorTickPosition-0.5*obj.MinorTickLength;
                R2 = obj.MinorTickPosition+0.5*obj.MinorTickLength;
            end
            
            % Creating tick end points
            P1 = [R1*cos(obj.ThtMinorTick); R1*sin(obj.ThtMinorTick)];
            P2 = [R2*cos(obj.ThtMinorTick); R2*sin(obj.ThtMinorTick)];
            
            % Creating auxiliary color
            if strcmpi(obj.Enable,'on')
                ColorAux = obj.MinorTickColor;
            else
                ColorAux = obj.OffMinorTickColor;
            end

            % Checking for change in number of ticks
            if strcmpi(Option,'new')
                
                % Deleting old ticks
                for i = 1:length(obj.MinorTick)
                    delete(obj.MinorTick{i})
                end
                obj.MinorTick = cell(1,length(obj.ThtMinorTick));
                
                % Creating new ticks
                for i = 1:length(obj.MinorTick)
                    obj.MinorTick{i} = plot(obj.Axis,[P1(1,i) P2(1,i)],[P1(2,i) P2(2,i)]);
                    set(obj.MinorTick{i},...
                        'Color',ColorAux,...
                        'LineWidth',obj.MinorTickLineWidth);
                end
                
            else
                
                % Updating existing ticks
                for i = 1:length(obj.MinorTick)
                    obj.MinorTick{i}.XData = [P1(1,i) P2(1,i)];
                    obj.MinorTick{i}.YData = [P1(2,i) P2(2,i)];
                end
                
            end
            
        end
        
        function drawTick(obj,varargin)
            
            % Checking for update option
            Option = 'new';
            if nargin == 2
                if strcmpi(varargin{1},'new')
                    Option = 'new';
                elseif strcmpi(varargin{1},'updated')
                    Option = 'updated';
                end
            end
            
            % Creating tick end point radius
            if strcmpi(obj.TickDir,'in')
                R1 = obj.TickPosition-obj.TickLength;
                R2 = obj.TickPosition;
            elseif strcmpi(obj.TickDir,'out')
                R1 = obj.TickPosition;
                R2 = obj.TickPosition+obj.TickLength;
            else
                R1 = obj.TickPosition-0.5*obj.TickLength;
                R2 = obj.TickPosition+0.5*obj.TickLength;
            end
            
            % Creating tick end points
            P1 = [R1*cos(obj.ThtTick); R1*sin(obj.ThtTick)];
            P2 = [R2*cos(obj.ThtTick); R2*sin(obj.ThtTick)];
            
            % Creating auxiliary color
            if strcmpi(obj.Enable,'on')
                ColorAux = obj.TickColor;
            else
                ColorAux = obj.OffTickColor;
            end

            % Checking for change in number of ticks
            if strcmpi(Option,'new')
                
                % Deleting old ticks
                for i = 1:length(obj.Tick)
                    delete(obj.Tick{i})
                end
                obj.Tick = cell(1,length(obj.ThtTick));
                
                % Creating new ticks
                for i = 1:length(obj.Tick)
                    obj.Tick{i} = plot(obj.Axis,[P1(1,i) P2(1,i)],[P1(2,i) P2(2,i)]);
                    set(obj.Tick{i},...
                        'Color',ColorAux,...
                        'LineWidth',obj.TickLineWidth);
                end
                
            else
                
                % Updating existing ticks
                for i = 1:length(obj.Tick)
                    obj.Tick{i}.XData = [P1(1,i) P2(1,i)];
                    obj.Tick{i}.YData = [P1(2,i) P2(2,i)];
                end
                
            end
            
        end
        
        function drawLabel(obj,varargin)
            
            % Checking for update option
            Option = 'new';
            if nargin == 2
                if strcmpi(varargin{1},'new')
                    Option = 'new';
                elseif strcmpi(varargin{1},'updated')
                    Option = 'updated';
                end
            end
            
            % Creating label radius
            if strcmpi(obj.LabelDir,'in')
                if strcmpi(obj.TickDir,'in')
                    R = obj.TickPosition-obj.TickLength-0.1;
                elseif strcmpi(obj.TickDir,'out')
                    R = obj.TickPosition-0.1;
                else
                    R = obj.TickPosition-0.5*obj.TickLength-0.1;
                end
            elseif strcmpi(obj.LabelDir,'out')
                d = 0.1 + 0.005*(obj.LabelFontSize-10);
                if strcmpi(obj.TickDir,'in')
                    R = obj.TickPosition+d;
                elseif strcmpi(obj.TickDir,'out')
                    R = obj.TickPosition+obj.TickLength+d;
                else
                    R = obj.TickPosition+0.5*obj.TickLength+d;
                end
            else
                R = obj.TickPosition;
            end
            
            % Creating label point
            P = [R*cos(obj.ThtTick); R*sin(obj.ThtTick)];

            % Creating auxiliary color
            if strcmpi(obj.Enable,'on')
                ColorAux = obj.LabelColor;
            else
                ColorAux = obj.OffLabelColor;
            end

            % Checking for change in number of labels
            if strcmpi(Option,'new')
                
                % Deleting old labels
                for i = 1:length(obj.Label)
                    delete(obj.Label{i})
                end
                obj.Label = cell(1,length(obj.ThtTick));
                
                % Updating label string
                StringAux = cell(1,length(obj.DialValues));
                for i = 1:length(obj.DialValues)
                    StringAux{i} = num2str(obj.DialValues(i));
                end
                obj.LabelString = StringAux;
                
                % Creating label string valules and new labels
                axes(obj.Axis)
                for i = 1:length(obj.Label)
                    obj.Label{i} = text(P(1,i),P(2,i),obj.LabelString{i});
                    set(obj.Label{i},...
                        'Color',ColorAux,...
                        'FontName',obj.LabelFontName,...
                        'FontSize',obj.LabelFontSize,...
                        'FontWeight',obj.LabelFontWeight,...
                        'HorizontalAlignment','center');
                end
                
            else
                
                % Updating existing labels
                for i = 1:length(obj.Label)
                    obj.Label{i}.Position = [P(1,i) P(2,i) 0];
                end
                
            end
                
        end
        
        function drawNeedle(obj)
            
            % Creating needle points
            if isempty(obj.NeedlePoints)
                if numel(obj.NeedleWidth) == 1
                    obj.NeedlePoints = [...
                        -obj.NeedleWidth/2            0         ; ...
                         obj.NeedleWidth/2            0         ; ...
                                 0            obj.NeedleLength ];
                else
                    obj.NeedlePoints = [...
                        -obj.NeedleWidth(1)/2         0         ; ...
                         obj.NeedleWidth(1)/2         0         ; ...
                         obj.NeedleWidth(2)/2 obj.NeedleLength  ;...
                        -obj.NeedleWidth(2)/2 obj.NeedleLength ];
                end
            end
            
            % Checking for needle
            if ~isempty(obj.Needle)
                delete(obj.Needle)
            end

            % Creating auxiliary color
            if strcmpi(obj.Enable,'on')
                ColorAux1 = obj.NeedleColor;
                ColorAux2 = obj.NeedleEdgeColor;
            else
                ColorAux1 = obj.OffNeedleColor;
                ColorAux2 = obj.OffNeedleEdgeColor;
            end

            % Creating needle
            axes(obj.Axis)
            obj.Needle = patch(...
                obj.NeedlePoints(:,1),...
                obj.NeedlePoints(:,2),...
                ColorAux1);
            set(obj.Needle,...
                'FaceLighting','flat',...
                'EdgeColor',ColorAux2,...
                'LineWidth',obj.NeedleLineWidth);
            
            % Adding needle object to transofmation group
            obj.NeedleGrp = hgtransform('Parent',obj.Axis);
            set(obj.Needle,'Parent',obj.NeedleGrp)
            
            % Moving needle
            obj.moveNeedle(obj.NeedlePosition)
            
        end
        
        function moveNeedle(obj,Value)
            
            % Calculating needle position
            ThtNeedle = interp1(obj.DialMinorValues,obj.ThtMinorTick,Value)-pi/2;
            
            % Defining rotation matrix
            T01 = [cos(ThtNeedle) -sin(ThtNeedle)        0             0 ;...
                   sin(ThtNeedle)  cos(ThtNeedle)        0             0 ;...
                          0               0              1             0 ;...
                          0               0              0             1 ];
                   
            % Defining translation matrix
            T12 = [    1         0         0         0          ;...
                       0         1         0   obj.NeedleOffset ;...
                       0         0         1         0          ;...
                       0         0         0         1          ];
                   
            % Transforming needle object
            obj.NeedleGrp.Matrix = T01*T12;
                   
        end
        
        function drawDisplay(obj,varargin)
            
            % Checking for update option
            Option = 'new';
            if nargin == 2
                if strcmpi(varargin{1},'new')
                    Option = 'new';
                elseif strcmpi(varargin{1},'updated')
                    Option = 'updated';
                elseif strcmpi(varargin{1},'value')
                    Option = 'value';
                end
            end
            
            % Assigning number formating
            NumFloat = max(0,3-length(num2str(floor(max(abs(obj.DialValues(end)))))));
            NumFormat = ['%0.' int2str(NumFloat) 'f'];
            
            % Assigning variable for conciseness
            PosAux = obj.DisplayPosition;
            
            % Calculating display center position
            PtDisp0 = PosAux(1)*[cosd(PosAux(2)-90); sind(PosAux(2)-90)];
                                         
            % Calculating display corners
            PtDisp = PtDisp0*ones(1,4) + ...
                0.5*[-PosAux(3) +PosAux(3) +PosAux(3) -PosAux(3); ...
                     -PosAux(4) -PosAux(4) +PosAux(4) +PosAux(4)];
            
            % Creating auxiliary color
            if strcmpi(obj.Enable,'on')
                ColorAux1 = obj.DisplayBackgroundColor;
                ColorAux2 = obj.DisplayEdgeColor;
                ColorAux3 = obj.DisplayForegroundColor;
            else
                ColorAux1 = obj.OffDisplayBackgroundColor;
                ColorAux2 = obj.OffDisplayEdgeColor;
                ColorAux3 = obj.OffDisplayForegroundColor;
            end

            % Checking for existing display
            switch Option
                case 'new'
                    axes(obj.Axis)
                    obj.Display{1} = patch(PtDisp(1,:),PtDisp(2,:),ColorAux1);
                    set(obj.Display{1},...
                        'FaceLighting','flat',...
                        'EdgeColor',ColorAux2,...
                        'LineWidth',obj.DisplayLineWidth);
                    obj.Display{2} = text(PtDisp0(1),PtDisp0(2),num2str(obj.NeedlePosition,NumFormat));
                    set(obj.Display{2},...
                        'Color',ColorAux3,...
                        'FontName',obj.DisplayFontName,...
                        'FontSize',obj.DisplayFontSize,...
                        'FontWeight',obj.DisplayFontWeight,...
                        'HorizontalAlignment','center');
                case 'updated'
                    obj.Display{1}.Vertices = PtDisp';
                    obj.Display{2}.Position = [PtDisp0(1) PtDisp0(2) 0];
                    obj.Display{2}.String = num2str(obj.NeedlePosition,NumFormat);
                case 'value'
                    obj.Display{2}.String = num2str(obj.NeedlePosition,NumFormat);
            end
            
        end
    
        function drawUnit(obj,varargin)
            
            % Checking for update option
            Option = 'new';
            if nargin == 2
                if strcmpi(varargin{1},'new')
                    Option = 'new';
                elseif strcmpi(varargin{1},'updated')
                    Option = 'updated';
                end
            end
            
            % Assigning variable for conciseness
            PosAux = obj.UnitPosition;
            
            % Calculating display center position
            PtUnit = PosAux(1)*[cosd(PosAux(2)-90); sind(PosAux(2)-90)];
                                         
            % Creating auxiliary color
            if strcmpi(obj.Enable,'on')
                ColorAux = obj.UnitColor;
            else
                ColorAux = obj.OffUnitColor;
            end

            % Checking for existing display
            if strcmpi(Option,'new')
                axes(obj.Axis)
                obj.Unit = text(PtUnit(1),PtUnit(2),num2str(obj.UnitString));
                set(obj.Unit,...
                    'Color',ColorAux,...
                    'FontName',obj.UnitFontName,...
                    'FontSize',obj.UnitFontSize,...
                    'FontWeight',obj.UnitFontWeight,...
                    'HorizontalAlignment','center');
            else
                obj.Unit.Position = [PtUnit(1) PtUnit(2) 0];
            end
            
        end
        
        
% =========================================================================           
%                  DEFINING OTHER USEFUL GAGE METHODS    
% =========================================================================           
    
        function deleteaxis(obj,~,~,~)
            
            % Deleting base axis
            obj.delete
        
        end
        
        function delete(obj)
            
            % Deleting base axis
            if ishandle(obj.Axis)
                delete(obj.Axis)
            end
        
        end
        
        function setdisp(~)
            
            disp('                    Parent: {}')
            disp('                       Box: {''on''  ''off''}')
            disp('                  Position: {}')
            disp('                     Units: {''pixels''  ''normalized''}')
            disp('                   Visible: {''on''  ''off''}')
            disp('                    Enable: {''on''  ''off''}')
            disp('                 BackColor: {}')
            disp('             BackEdgeColor: {}')
            disp('             BackLineWidth: {}')
            disp('               BackVisible: {''on''  ''off''}')
            disp('                BandValues: {}')
            disp('                 BandColor: {}')
            disp('              BandPosition: [InnerRadius  OuterRadius]')
            disp('               BandVisible: {''on''  ''off''}')
            disp('                 DialAngle: {}')
            disp('                DialValues: {}')
            disp('           DialMinorValues: {}')
            disp('                 DialColor: {}')
            disp('             DialLineWidth: {}')
            disp('               DialVisible: {''on''  ''off''}')
            disp('               LabelString: {}')
            disp('             LabelFontName: {}')
            disp('             LabelFontSize: {}')
            disp('           LabelFontWeight: {''normal''  ''bold''}')
            disp('                LabelColor: {}')
            disp('                  LabelDir: {''in''  ''center''  ''out''}')
            disp('          LabelOrientation: {''horizontal''  ''radial''  ''tangential''}')
            disp('              LabelVisible: {''on''  ''off''}')
            disp('                 TickColor: {}')
            disp('                   TickDir: {''in''  ''center''  ''out''}')
            disp('                TickLength: {}')
            disp('             TickLineWidth: {}')
            disp('              TickPosition: [Radius]')
            disp('               TickVisible: {''on''  ''off''}')
            disp('            MinorTickColor: {}')
            disp('              MinorTickDir: {''in''  ''center''  ''out''}')
            disp('           MinorTickLength: {}')
            disp('        MinorTickLineWidth: {}')
            disp('         MinorTickPosition: [Radius]')
            disp('          MinorTickVisible: {''on''  ''off''}')
            disp('                  KnobSize: {}')
            disp('                 KnobColor: {}')
            disp('             KnobEdgeColor: {}')
            disp('             KnobLineWidth: {}')
            disp('               NeedleColor: {}')
            disp('           NeedleEdgeColor: {}')
            disp('           NeedleLineWidth: {}')
            disp('              NeedleLength: {}')
            disp('               NeedleWidth: {}')
            disp('              NeedleOffset: {}')
            disp('              NeedlePoints: {}')
            disp('            NeedlePosition: [Value]')
            disp('           DisplayFontName: {}')
            disp('           DisplayFontSize: {}')
            disp('         DisplayFontWeight: {''normal''  ''bold''}')
            disp('    DisplayForegroundColor: {}')
            disp('    DisplayBackgroundColor: {}')
            disp('          DisplayEdgeColor: {}')
            disp('          DisplayLineWidth: {}')
            disp('           DisplayPosition: [Radius Angle BoxWidth BoxHeight]')
            disp('            DisplayVisible: {}')
            disp('                UnitString: {}')
            disp('              UnitFontName: {}')
            disp('              UnitFontSize: {}')
            disp('            UnitFontWeight: {''normal''  ''bold''}')
            disp('                 UnitColor: {}')
            disp('              UnitPosition: [Radius Angle]')
            disp('               UnitVisible: {''on''  ''off''}')
            
        end

        function putProperties(obj,nProp,PropList,iStart)
            
            % Looping through property list
            for i = iStart:2:nProp-1
                switch lower(PropList{i})
                    case 'parent'
                        % Taken care of in the constructor method
                    case 'box'
                        obj.Box = PropList{i+1};
                    case 'position'
                        % Taken care of in the constructor method
                    case 'units'
                        % Taken care of in the constructor method
                    case 'visible'
                        obj.Visible = PropList{i+1};
                    case 'enable'
                        % Taken care of in the constructor method
                    case 'backcolor'
                        obj.BackColor = PropList{i+1};
                    case 'backedgecolor'
                        obj.BackEdgeColor = PropList{i+1};
                    case 'backlinewidth'
                        obj.BackLineWidth = PropList{i+1};
                    case 'backvisible'
                        obj.BackVisible = PropList{i+1};
                    case 'bandvalues'
                        obj.BandValues = PropList{i+1};
                    case 'bandcolor'
                        obj.BandColor = PropList{i+1};
                    case 'bandposition'
                        obj.BandPosition = PropList{i+1};
                    case 'bandvisible'
                        obj.BandVisible = PropList{i+1};
                    case 'dialangle'
                        obj.DialAngle = PropList{i+1};
                    case 'dialvalues'
                        obj.DialValues = PropList{i+1};
                    case 'dialminorvalues'
                        obj.DialMinorValues = PropList{i+1};
                    case 'dialcolor'
                        obj.DialColor = PropList{i+1};
                    case 'diallinewidth'
                        obj.DialLineWidth = PropList{i+1};
                    case 'dialvisible'
                        obj.DialVisible = PropList{i+1};
                    case 'labelstring'
                        % Taken care of in the constructor method
                    case 'labelfontname'
                        obj.LabelFontName = PropList{i+1};
                    case 'labelfontsize'
                        obj.LabelFontSize = PropList{i+1};
                    case 'labelfontweight'
                        obj.LabelFontWeight = PropList{i+1};
                    case 'labelcolor'
                        obj.LabelColor = PropList{i+1};
                    case 'labeldir'
                        obj.LabelDir = PropList{i+1};
                    case 'labelorientation'
                        obj.LabelOrientation = PropList{i+1};
                    case 'labelvisible'
                        obj.LabelVisible = PropList{i+1};
                    case 'tickcolor'
                        obj.TickColor = PropList{i+1};
                    case 'tickdir'
                        obj.TickDir = PropList{i+1};
                    case 'ticklength'
                        obj.TickLength = PropList{i+1};
                    case 'ticklinewidth'
                        obj.TickLineWidth = PropList{i+1};
                    case 'tickposition'
                        obj.TickPosition = PropList{i+1};
                    case 'tickvisible'
                        obj.TickVisible = PropList{i+1};
                    case 'minortickcolor'
                        obj.MinorTickColor = PropList{i+1};
                    case 'minortickdir'
                        obj.MinorTickDir = PropList{i+1};
                    case 'minorticklength'
                        obj.MinorTickLength = PropList{i+1};
                    case 'minorticklinewidth'
                        obj.MinorTickLineWidth = PropList{i+1};
                    case 'minortickposition'
                        obj.MinorTickPosition = PropList{i+1};
                    case 'minortickvisible'
                        obj.MinorTickVisible = PropList{i+1};
                    case 'knobsize'
                        obj.KnobSize = PropList{i+1};
                    case 'knobcolor'
                        obj.KnobColor = PropList{i+1};
                    case 'knobedgecolor'
                        obj.KnobEdgeColor = PropList{i+1};
                    case 'knoblinewidth'
                        obj.KnobLineWidth = PropList{i+1};
                    case 'needlecolor'
                        obj.NeedleColor = PropList{i+1};
                    case 'needleedgecolor'
                        obj.NeedleEdgeColor = PropList{i+1};
                    case 'needlelinewidth'
                        obj.NeedleLineWidth = PropList{i+1};
                    case 'needlelength'
                        obj.NeedleLength = PropList{i+1};
                    case 'needlewidth'
                        obj.NeedleWidth = PropList{i+1};
                    case 'needleoffset'
                        obj.NeedleOffset = PropList{i+1};
                    case 'needlepoints'
                        obj.NeedlePoints = PropList{i+1};
                    case 'needleposition'
                        obj.NeedlePosition = PropList{i+1};
                    case 'displayfontname'
                        obj.DisplayFontName = PropList{i+1};
                    case 'displayfontsize'
                        obj.DisplayFontSize = PropList{i+1};
                    case 'displayfontweight'
                        obj.DisplayFontWeight = PropList{i+1};
                    case 'displayforegroundcolor'
                        obj.DisplayForegroundColor = PropList{i+1};
                    case 'displaybackgroundcolor'
                        obj.DisplayBackgroundColor = PropList{i+1};
                    case 'displayedgecolor'
                        obj.DisplayEdgeColor = PropList{i+1};
                    case 'displaylinewidth'
                        obj.DisplayLineWidth = PropList{i+1};
                    case 'displayposition'
                        obj.DisplayPosition = PropList{i+1};
                    case 'displayvisible'
                        obj.DisplayVisible = PropList{i+1};
                    case 'unitstring'
                        obj.UnitString = PropList{i+1};
                    case 'unitfontname'
                        obj.UnitFontName = PropList{i+1};
                    case 'unitfontsize'
                        obj.UnitFontSize = PropList{i+1};
                    case 'unitfontweight'
                        obj.UnitFontWeight = PropList{i+1};
                    case 'unitcolor'
                        obj.UnitColor = PropList{i+1};
                    case 'unitposition'
                        obj.UnitPosition = PropList{i+1};
                    case 'unitvisible'
                        obj.UnitVisible = PropList{i+1};
                    otherwise
                        error(['Unknown property : ' PropList{i}]);
                end
            end
            
        end
        
    end
    
    
% ========================================================================= 
%                      DEFINING AUXILIARY FUNCTIONS
% =========================================================================           
    
    methods (Static = true, Hidden = true)
        
        function [AngleArray,BandArray,ValueArray,MinorValuesArray] = ...
                calculateAngles(Angle,BandValues,Values,MinorValues)
            
            % Calculating dial start and end angles
            ThtStart = -pi/180*(-Angle(2)+90);
            ThtEnd   =  pi/180*(270+Angle(1));

            % Creating dial (unit radius)
            AngleArray = ThtStart:pi/180:ThtEnd;
            
            % Creating color band
            BandArray = cell(1,length(BandValues)-1);
            for i = 1:length(BandArray)
                ThtBandStart = ThtEnd + (ThtStart-ThtEnd)/(Values(end)-Values(1))*(BandValues(i)-Values(1));
                ThtBandEnd = ThtEnd + (ThtStart-ThtEnd)/(Values(end)-Values(1))*(BandValues(i+1)-Values(1));
                nBand = ceil((ThtBandStart-ThtBandEnd)/(pi/180));
                BandArray{i} = ThtBandStart:-(ThtBandStart-ThtBandEnd)/(nBand-1):ThtBandEnd;
            end
            % Creating angle positions for major ticks and labels
            ValueArray = ThtEnd + (ThtStart-ThtEnd)/(Values(end)-Values(1))*(Values-Values(1));

            % Creating angle positions for minor ticks
            MinorValuesArray = ThtEnd + (ThtStart-ThtEnd)/(Values(end)-Values(1))*(MinorValues-Values(1));

        end
        
    end
    
end