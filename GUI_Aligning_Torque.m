% GUI ALIGNING TORQUE
% Written by: JoshTheEngineer
% Started: 09/06/17
% Updated: 09/06/17 - Started code
%                   - Works as intended with all plotting capabilities

function varargout = GUI_Aligning_Torque(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Aligning_Torque_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Aligning_Torque_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% --- Executes just before GUI_Aligning_Torque is made visible.
function GUI_Aligning_Torque_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = GUI_Aligning_Torque_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% INITIALIZATION

% Call the SOLVE function
SOLVE(handles);

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% --------------------------- INITIALIZATION ---------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

% EDIT ----------------------- Slip Angle ---------------------------------
function editSlipAngle_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% LIST ----------------------- Select Plot --------------------------------
function listPlot_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------- CALLBACKS ------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

% EDIT ----------------------- Slip Angle ---------------------------------
function editSlipAngle_Callback(hObject, eventdata, handles)
% =========================================================================
% - Read in the slip angle from the edit text box
% - Compute the solution
% - Plot based on the selection in the list box
% =========================================================================

% Call the SOLVE function
SOLVE(handles);

% LIST ----------------------- Select Plot --------------------------------
function listPlot_Callback(hObject, eventdata, handles)
% =========================================================================
% - Call the PLOT function
% =========================================================================

% Call the PLOT function
PLOT(handles);

% FUNCTION ------------------- S O L V E ----------------------------------
function [] = SOLVE(handles)
% =========================================================================
% - Read in the slip angle from the edit text box
% - Compute the solution
% - Plot based on the selection in the list box
% =========================================================================

% Define the slip angle
SAval = str2double(get(handles.editSlipAngle,'String'));                    % Read in user defined value
SA    = (0:0.5:90)';                                                        % Define total array over valid range

% Line-circle intersection
r  = 1;                                                                     % Radius of lateral force potential
b  = tand(SA);                                                              % Y-intercept of triangle hypotenuse line
m  = b;                                                                     % Slope of triangle hypotenuse line
A  = 1 + b.^2;                                                              % A component of quadratic equation
B  = 2.*m.*b;                                                               % B component of quadratic equation
C  = b.^2 - r.^2;                                                           % C component of quadratic equation

Xc = (-B + sqrt(B.^2-(4.*A.*C)))./(2.*A);                                   % Solve quadratic equation for Xc
Yc = m.*Xc + b;                                                             % Linear equation to solve for Yc

% Triangle centroid, lateral force, and moment
Xbar_T = ((2/3).*(Xc+1))-1;                                                 % Centroid of triangle part of footprint
Fy_T   = (1/2).*(1+Xc).*Yc;                                                 % Triangle total lateral force [N] = 0.5*base*height
Mz_T   = Fy_T.*Xbar_T;                                                      % Triangle aligning moment [N*m]

% Circular centroid, lateral force, and moment
num    = ((1-Xc.^2).^(3/2))./3;                                             % Integral of x*sqrt(1-x^2)dx from [Xc,1]
den    = asin(1)/2 - asin(Xc)./2 - (1/2).*Xc.*((1-Xc.^2).^(1/2));           % Integral of sqrt(1-x^2)dx from [Xc,1]
Xbar_C = num./den;                                                          % Centroid of a circle segment
Fy_C   = den;                                                               % Circle total lateral force [N]
Mz_C   = Fy_C.*Xbar_C;                                                      % Circle aligning moment [N*m]

% Total lateral force, pneumatic trail, and aligning moment
Fy = Fy_T + Fy_C;                                                           % Total lateral force from footprint [N]
Mz = Mz_T + Mz_C;                                                           % Total aligning moment from footprint [N*m]
t  = Mz./Fy;                                                                % Pneumatic trail [m]

% Find the values corresponding to the desired SAval
dSA   = SA(2)-SA(1);                                                        % Slip angle step in array
indSA = find(SA > SAval-dSA/2 & SA < SAval+dSA/2);                          % Index of selected slip angle

% Save the solution data to an array
solArray = [SA ...                                                          % Construct the solution array
            Xc     Yc ...
            Xbar_T Xbar_C ...
            Fy_T   Fy_C ...
            Mz_T   Mz_C ...
            Fy Mz t];
solArray(1,5)  = 1;
solArray(1,9)  = 0;
solArray(1,11) = 0;
solArray(1,12) = 1/3;
assignin('base','solArray',solArray);
assignin('base','indSA',indSA);

% Call the PLOT function
PLOT(handles);

% FUNCTION -------------------- P L O T -----------------------------------
function [] = PLOT(handles)
% =========================================================================
% - Read in the solution array
% - Get list box selection for what to plot
% - Plot the data
% =========================================================================

% Get the list box selection
listVal = get(handles.listPlot,'Value');

% Read in the solution array
solArray = evalin('base','solArray');                                       % Solution array
indSA    = evalin('base','indSA');                                          % Index corresponding to the user-defined SA

% Plot the friction circle
axes(handles.plotFrictionCircle);
cla reset;
hold on; grid on;                                                           % Get ready for plotting
theta = linspace(0,pi,100)';                                                % Define theta angle for plotting lateral force potential
plot(cos(theta),sin(theta),'k--','LineWidth',2);                            % Plot the lateral force potential circle
xlim([-1 1]);                                                               % X-axis limits
plot(solArray(indSA,2),solArray(indSA,3),'ro',...                           % Plot the line-circle intersection point
     'MarkerFaceColor','k','MarkerEdgeColor','k');
plot([-1 solArray(indSA,2)],[0 solArray(indSA,3)],'b-','LineWidth',2);
XT_Fill = [-1 solArray(indSA,2) solArray(indSA,2) -1];
YT_Fill = [0  solArray(indSA,3)  0  0];
fill(XT_Fill,YT_Fill,'b');
plot([solArray(indSA,4) solArray(indSA,4)],[0 1],'k--','LineWidth',2);
XC_Fill = linspace(0,1-solArray(indSA,2))';
YC_Fill = sqrt(1-(XC_Fill+solArray(indSA,2)).^2);
XC_Fill = [XC_Fill; 0; 0];
YC_Fill = [YC_Fill; 0; solArray(indSA,3)];
fill(solArray(indSA,2)+XC_Fill,YC_Fill,'r');
plot([solArray(indSA,5) solArray(indSA,5)],[0 1],'k--','LineWidth',2);
plot([solArray(indSA,12) solArray(indSA,12)],[0 1],'k','LineWidth',4');
title('Tire Footprint');
xlabel('X [arb]');
ylabel('Y [arb]');
axis('equal');

% Plot based off of list selection
axes(handles.plotData);                                                     % Select the appropriate plot
if (listVal == 1)                                                           % :: Centroid vs. Alpha
    cla reset;                                                              % Clear and reset axes
    hold on; grid on;                                                       % Get ready for plotting
    p1 = plot(solArray(:,1),solArray(:,4),'b-','LineWidth',2);              % Plot the triangle centroid
    p2 = plot(solArray(:,1),solArray(:,5),'r-','LineWidth',2);              % Plot the circle centroid
    p3 = plot(solArray(indSA,1),solArray(indSA,4),'bo');                    % Plot the triangle point corresponding to SA
    p4 = plot(solArray(indSA,1),solArray(indSA,5),'ro');                    % Plot the circle point corresponding to SA
    set(p3,'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',10);    % Adjust the point's properties
    set(p4,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10);    % Adjust the point's properties
    title('Centroid Location vs. Slip Angle');                              % Title
    xlabel('Slip Angle [deg]');                                             % X-axis label
    ylabel('Centroid X Location [m]');                                      % Y-axis label
    legend([p1 p2],'Triangle','Circle');                                    % Legend
    ylim('auto');                                                           % Auto-scale the Y-axis
    xlim([0 90]);                                                           % Set the X-axis limits
elseif (listVal == 2)                                                       % :: Mz vs. Alpha (Individual)
    cla reset;                                                              % Clear and reset axes
    hold on; grid on;                                                       % Get ready for plotting
    p1 = plot(solArray(:,1),solArray(:,8),'b-','LineWidth',2);              % Plot the triangle Mz
    p2 = plot(solArray(:,1),solArray(:,9),'r-','LineWidth',2);              % Plot the circle Mz
    p3 = plot(solArray(indSA,1),solArray(indSA,8),'bo');                    % Plot the triangle point corresponding to SA
    p4 = plot(solArray(indSA,1),solArray(indSA,9),'ro');                    % Plot the circle point corresponding to SA
    set(p3,'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',10);    % Adjust the point's properties
    set(p4,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10);    % Adjust the point's properties
    title('Aligning Moment vs. Slip Angle');                                % Title
    xlabel('Slip Angle [deg]');                                             % X-axis label
    ylabel('Aligning Moment [N*m]');                                        % Y-axis label
    legend([p1 p2],'Triangle','Circle');                                    % Legend
    ylim('auto');                                                           % Auto-scale the Y-axis
    xlim([0 90]);                                                           % Set the X-axis limits
elseif (listVal == 3)                                                       % :: Fy and t vs. Alpha
    cla reset;                                                              % Clear and reset axes
    hold on; grid on;                                                       % Get ready for plotting
    [AX,H1,H2] = plotyy(solArray(:,1),solArray(:,10),...                    % Plot the Fy and t for both left and right Y-axes
                        solArray(:,1),solArray(:,12));
    set(H1,'Color','k','LineStyle','-','LineWidth',3);                      % Set the Fy plot color
    set(H2,'Color','k','LineStyle','-.','LineWidth',3);                     % Set the t plot color
    H3 = plot(solArray(:,1),solArray(:,11)*7,'b--','LineWidth',2);          % Plot the arb total Mz
    set(AX(1),'YColor','k');                                                % Set the left Y-axis color
    set(AX(2),'YColor','k');                                                % Set the right Y-axis color
    xlabel('Slip Angle [deg]');                                             % X-axis label
    ylabel(AX(1),'Lateral Force [N]','Color','k');                          % Left Y-axis label
    ylabel(AX(2),'Pneumatic Trail [m]','Color','k');                        % Right Y-axis label
    legend([H1 H2 H3],'Lateral Force (L)','Pneumatic Trail (R)',...         % Legend
                      'Arb Aligning Torque');
    xlim(AX(1),[0 90]);                                                     % Set X-axis limits for left plot
    xlim(AX(2),[0 90]);                                                     % Set X-axis limits for right plot
elseif (listVal == 4)                                                       % :: Mz vs. Alpha (One-Sided)
    cla reset;                                                              % Clear and reset axes
    hold on; grid on;                                                       % Get ready for plotting
    plot(solArray(:,1),solArray(:,11),'k-','LineWidth',3);                  % Plot the total aligning moment
    title('Aligning Moment vs. Slip Angle');                                % Title
    xlabel('Slip Angle [deg]');                                             % X-axis label
    ylabel('Aligning Torque [N*m]');                                        % Y-axis label
    ylim('auto');                                                           % Y-axis limits
    xlim([0 90]);                                                           % X-axis limits
elseif (listVal == 5)                                                       % :: Mz vs. Alpha (Two-Sided)
    cla reset;                                                              % Clear and reset axes
    hold on; grid on;                                                       % Get ready for plotting
    plot(solArray(:,1),solArray(:,11),'k-','LineWidth',3);                  % Plot the left-hand turn aligning moment
    plot(-solArray(:,1),-solArray(:,11),'k-','LineWidth',3);                % Plot the right-hand turn aligning moment
    xlabel('Slip Angle [deg]');                                             % X-axis label
    ylabel('Aligning Torque [N*m]');                                        % Y-axis label
    ylim('auto');                                                           % Y-axis limits
    xlim([-90 90]);                                                         % X-axis limits
end

% PUSH -------------------------- Exit ------------------------------------
function pushExit_Callback(hObject, eventdata, handles)
delete(handles.figureGUIAligningTorque);
clc;
