function varargout = plentyOfRope(varargin)
% PLENTYOFROPE MATLAB code for plentyOfRope.fig
%      PLENTYOFROPE, by itself, creates a new PLENTYOFROPE or raises the existing
%      singleton*.
%
%      H = PLENTYOFROPE returns the handle to a new PLENTYOFROPE or the handle to
%      the existing singleton*.
%
%      PLENTYOFROPE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLENTYOFROPE.M with the given input arguments.
%
%      PLENTYOFROPE('Property','Value',...) creates a new PLENTYOFROPE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plentyOfRope_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plentyOfRope_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plentyOfRope

% Last Modified by GUIDE v2.5 13-Jul-2015 04:30:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plentyOfRope_OpeningFcn, ...
                   'gui_OutputFcn',  @plentyOfRope_OutputFcn, ...
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
% End initialization code - DO NOT EDIT
end


% --- Executes just before plentyOfRope is made visible.
function plentyOfRope_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plentyOfRope (see VARARGIN)

% Choose default command line output for plentyOfRope
handles.output = hObject;

% UIWAIT makes plentyOfRope wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Initialize dynamics model parameters:

handles.tau = 100;
set(handles.input_tau,'String',handles.tau);
handles.difu = 30;
set(handles.input_difu,'String',handles.difu);
handles.temp = 293;
set(handles.input_temp,'String',handles.temp);
q = 1;
kB = 8.617e-5; % eV/K
handles.mu = handles.difu*q/kB/handles.temp;
set(handles.input_mu,'String',handles.mu);
handles.SRV = 1e6;
set(handles.input_SRV,'String',handles.SRV);
set(handles.menu_solver, 'Value', 1); % start with analytical model
solverOptions = get(handles.menu_solver, 'String');
handles.solver = solverOptions{1};
set(handles.menu_tauModel, 'Value', 1);
tauModelOptions = get(handles.menu_tauModel, 'String');
handles.tauModel = tauModelOptions{1};
set(handles.menu_difuModel, 'Value', 1);
difuModelOptions = get(handles.menu_difuModel, 'String');
handles.difuModel = difuModelOptions{1};
handles.effMassRatio    = 2;
set(handles.input_effMassRatio, 'String', handles.effMassRatio);
handles.majConc         = 1e16;
set(handles.input_majConc, 'String', handles.majConc);
handles.alpha = 292.4;
set(handles.input_alpha,'String',handles.alpha);
handles.R = 0.3;
set(handles.input_R,'String',handles.R);
handles.thick = 400;
set(handles.input_d,'String',handles.thick);

% Initialize pump model parameters:

set(handles.menu_pumpProfile,'Value',1);
pumpProfileOptions = get(handles.menu_pumpProfile, 'String');
handles.pumpProfile = pumpProfileOptions{1};
handles.N = 1;
set(handles.input_N,'String',handles.N);
handles.T = 1e-3;
set(handles.input_T,'String',handles.T);
handles.sigma = 1e-3;
set(handles.input_sigma,'String',handles.sigma);
handles.h = 0;
set(handles.input_h, 'String', handles.h);

% Initialize the data, offsets, and scaling:

handles.dataX           = NaN;
handles.dataY           = NaN;
handles.dataXShift      = 0;
set(handles.input_dataXShift, 'String', handles.dataXShift);
handles.dataYShift      = 0;
set(handles.input_dataYShift, 'String', handles.dataYShift);
handles.nomMajMob       = 10;
set(handles.input_nomMajMob, 'String', handles.nomMajMob);
set(handles.text_dataError, 'visible', 'off');
set(handles.input_datasetName, 'String', 'noname');
set(handles.indicate_fitNum, 'String', num2str(0,'%u'));

% Initialize the radio buttons:

set(handles.input_fit_tau,          'Value', 1);
set(handles.input_fit_difu,           'Value', 1);
set(handles.input_fit_SRV,            'Value', 1);
set(handles.input_fit_alpha,        'Value', 1);
set(handles.input_fit_N,            'Value', 1);
set(handles.input_fit_sigma,        'Value', 1);
set(handles.input_fit_h,            'Value', 1);
set(handles.input_fit_dataXShift,   'Value', 1);
set(handles.input_fit_dataYShift,   'Value', 1);
set(handles.input_scaleByMobility,  'Value', 0);
set(handles.input_plotProfile,      'Value', 0);
set(handles.input_plotCovariance,   'Value', 0);
set(handles.input_latchModel,       'Value', 0);

% Make sure sliders and text inputs are equal:

SRVMin = floor(log10(handles.SRV));
SRVMax = SRVMin + 1;
set(handles.slider_SRV,'Min',10^SRVMin,'Max',10^SRVMax,'Value',handles.SRV);
SRVMinLabel = sprintf('10^%d',SRVMin);
SRVMaxLabel = sprintf('10^%d',SRVMax);
set(handles.slider_SRV_min_label,'String',SRVMinLabel);
set(handles.slider_SRV_max_label,'String',SRVMaxLabel);
tauMin = floor(log10(handles.tau));
tauMax = tauMin + 1;
set(handles.slider_tau,'Min',10^tauMin,'Max',10^tauMax,'Value',handles.tau);
tauMinLabel = sprintf('10^%d',tauMin);
tauMaxLabel = sprintf('10^%d',tauMax);
set(handles.slider_tau_min_label,'String',tauMinLabel);
set(handles.slider_tau_max_label,'String',tauMaxLabel);
difuMin = floor(log10(handles.difu));
difuMax = difuMin + 1;
set(handles.slider_difu,'Min',10^difuMin,'Max',10^difuMax,'Value',handles.difu);
difuMinLabel = sprintf('10^%.0f',difuMin);
difuMaxLabel = sprintf('10^%.0f',difuMax);
set(handles.slider_difu_min_label,'String',difuMinLabel);
set(handles.slider_difu_max_label,'String',difuMaxLabel);

% Initialize the main axes:

axes(handles.mainAxes)
xLim = [0 10];
yLim = [1 10];
set(handles.mainAxes,'XLim',xLim,'YLim',yLim);
set(handles.input_xAxesMin,'String',xLim(1));
set(handles.input_xAxesMax,'String',xLim(2));
set(handles.input_yAxesMin,'String',yLim(1));
set(handles.input_yAxesMax,'String',yLim(2));
c_ = get(handles.mainAxes,'children');
delete(c_)
hold all
set(handles.mainAxes,'YScale','log');
set(handles.input_logScale,'Value',1);
set(handles.mainAxes,'yLimMode','auto');
set(handles.input_yAutoScale,'Value',1);
set(handles.mainAxes,'xLimMode','auto');
set(handles.input_xAutoScale,'Value',1);

% Initialize the colored boxes:

set(handles.colorBox1, 'box', 'on', 'ytick', [], 'xtick', []);
set(handles.colorBox2, 'box', 'on', 'ytick', [], 'xtick', []);
set(handles.colorBox3, 'box', 'on', 'ytick', [], 'xtick', []);

% Update handles structure
guidata(hObject, handles);
end

% --- Outputs from this function are returned to the command line.
function varargout = plentyOfRope_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Updates all the input fields in the GUI ---------------------------
function updateInputs(handles)

% Update tau input and slider:

set(handles.input_tau,'String',sprintf('%.3e',handles.tau));
tauMin = floor(log10(handles.tau));
tauMax = tauMin + 1;
set(handles.slider_tau,'Min',10^tauMin,'Max',10^tauMax,'Value',handles.tau);
tauMinLabel = sprintf('10^%d',tauMin);
tauMaxLabel = sprintf('10^%d',tauMax);
set(handles.slider_tau_min_label,'String',tauMinLabel);
set(handles.slider_tau_max_label,'String',tauMaxLabel);

% Update diffusivity input and slider:

set(handles.input_difu,'String',sprintf('%.3f',handles.difu));
difuMin = floor(log10(handles.difu));
difuMax = difuMin + 1;
set(handles.slider_difu,'Min',10^difuMin,'Max',10^difuMax,'Value',handles.difu);
difuMinLabel = sprintf('10^%.0f',difuMin);
difuMaxLabel = sprintf('10^%.0f',difuMax);
set(handles.slider_difu_min_label,'String',difuMinLabel);
set(handles.slider_difu_max_label,'String',difuMaxLabel);

% Update SRV input and slider:

set(handles.input_SRV,'String',sprintf('%.3e',handles.SRV));
SRVMin = floor(log10(handles.SRV));
SRVMax = SRVMin + 1;
set(handles.slider_SRV,'Min',10^SRVMin,'Max',10^SRVMax,'Value',handles.SRV);
SRVMinLabel = sprintf('10^%d',SRVMin);
SRVMaxLabel = sprintf('10^%d',SRVMax);
set(handles.slider_SRV_min_label,'String',SRVMinLabel);
set(handles.slider_SRV_max_label,'String',SRVMaxLabel);

% Update other parameters:

set(handles.input_mu,           'String', sprintf('%.1f', handles.mu));
set(handles.input_temp,         'String', sprintf('%.1f', handles.temp));
set(handles.input_effMassRatio, 'String', sprintf('%.1f', handles.effMassRatio));
set(handles.input_majConc,      'String', sprintf('%.2e', handles.majConc));
set(handles.input_alpha,        'String', sprintf('%.3e', handles.alpha));
set(handles.input_R,            'String', handles.R);
set(handles.input_d,            'String', handles.thick);
set(handles.input_N,            'String', sprintf('%.2e', handles.N));
set(handles.input_T,            'String', handles.T);
set(handles.input_sigma,        'String', sprintf('%.2e', handles.sigma));
set(handles.input_h,            'String', sprintf('%.2e', handles.h));
set(handles.input_dataXShift,   'String', sprintf('%.2e', handles.dataXShift));
set(handles.input_dataYShift,   'String', sprintf('%.2e', handles.dataYShift));
set(handles.input_nomMajMob,    'String', sprintf('%.2e', handles.nomMajMob));

% Update injection dependence models

tauModelOptions = get(handles.menu_tauModel, 'String');
tf = strcmp(handles.tauModel, tauModelOptions);
tf = tf'*(1:length(tauModelOptions))';
set(handles.menu_tauModel, 'Value', tf);

difuModelOptions = get(handles.menu_difuModel, 'String');
tf = strcmp(handles.difuModel, difuModelOptions);
tf = tf'*(1:length(difuModelOptions))';
set(handles.menu_difuModel, 'Value', tf);

% Update plot axes control:

xLim = get(handles.mainAxes,'XLim');
set(handles.input_xAxesMin,'String',xLim(1));
set(handles.input_xAxesMax,'String',xLim(2));
yLim = get(handles.mainAxes,'YLim');
set(handles.input_yAxesMin,'String',yLim(1));
set(handles.input_yAxesMax,'String',yLim(2));

end


% --- Plots the model function using the current parameters -------------
function plotModel(handles)
% Evaluate the model using the current parameters, and over a range
% determined by the axes limits. Also calls plotData, in case the data is
% affected by the model (as in the case of mobility scaling)

if (get(handles.input_latchModel, 'Value'))
    return
end

plotData(handles);

axes(handles.mainAxes);

xLim = get(handles.mainAxes, 'XLim');
t = linspace(xLim(1), xLim(2), 1000); % timeseries 
tDiscrete = linspace(xLim(1), xLim(2), 5); 
% discrete points in time to plot the spatial profiles
y = linspace(-handles.thick/2, handles.thick/2, 100); 
% points in space through depth of wafer

P = [handles.tau, handles.SRV, handles.thick, handles.alpha, handles.R,...
    handles.difu, handles.N, handles.sigma, handles.T];

injDep = struct;
% structure to pass injection dependence info to solveFC
injDep.tauModel     = handles.tauModel;
injDep.difuModel    = handles.difuModel;
injDep.effMassRatio = handles.effMassRatio;
injDep.majConc      = handles.majConc;
injDep.temp         = handles.temp;

wait_ = [];
if strcmp(handles.solver, 'Numerical')
    f_ = get(handles.mainAxes, 'parent');
    wait_ = annotation(f_, 'textbox', [.15 .5 .2 .1], 'String',...
        'Calculating ...', 'EdgeColor', 'r', 'Color', 'r', 'BackgroundColor', 'w',...
        'FontSize', 25);
    pause(.05)
end
% throw up a "wait" signal for numerical solver

[nAvg, nProfile, mode_info] = solveFC(t, tDiscrete, y, handles.pumpProfile, ...
    P, handles.h, handles.solver, injDep, true);
% call solveFC

if ~isempty(wait_)
    delete(wait_)
end

% Find the plotted model lineseries and replace it:
axes(handles.mainAxes);
hM_ = findobj('tag', 'modelPlotLine'); 
hMP_ = findobj('tag', 'prevModelPlotLine');
if ~isempty(hMP_)
    delete(hMP_)
end
if ~isempty(hM_)
    set(hM_, 'linestyle', '--', 'linewidth', 1, 'tag', 'prevModelPlotLine')
end
plot(t, nAvg, 'k-', 'tag', 'modelPlotLine', 'linewidth', 2);

xlabel('time (ns)')
ylabel('average carrier concentraion n (cm^{-3})')

grid on
box on

% Plot excess carrier depth distribution:

fDepthDistrib_ = findobj('tag', 'fDepthDistrib');
if get(handles.input_plotProfile, 'Value')
    if isempty(fDepthDistrib_)
        fDepthDistrib_ = figure('name', 'Model n(x,t)', 'tag', 'fDepthDistrib',...
            'Position', [800 70 420 350]);
        a_ = axes('parent', fDepthDistrib_); 
        axes(a_); hold all
    else
        delete(get(fDepthDistrib_, 'children'));
        a_ = axes('parent', fDepthDistrib_);
        axes(a_); hold all
    end
    for i = 1:length(tDiscrete)
        plot(y, nProfile(:,i), 'displayname', sprintf('t=%.2e \ns', tDiscrete(i)),...
            'parent', a_);
    end
    legend(a_, 'show')
    xlabel('Position (\mum)')
    ylabel('local carrier concentration n (cm^{-3})')
else
    delete(fDepthDistrib_);
end

% Update the data table:

fDataTables_ = findobj('tag', 'fDataTables');
if isempty(fDataTables_)
    [tableModel_, tableData_] = makeTableFig;
else
    tableModel_ = findobj('tag', 'tableModel');
end
set(tableModel_, 'Data', [t' nAvg']);

% Update the decay mode plot, if analytical solver was used:

fDisplayModes_ = findobj('tag', 'fDisplayModes');
if ~isempty(fDisplayModes_)
    delete(fDisplayModes_);
end
% delete pre-existing figure
if strcmp(handles.solver, 'Analytical')
    fDisplayModes_ = figure('tag', 'fDisplayModes', 'name', 'analytical model decay modes',...
        'position', [430 600 470 380]);
    tauAlphaN = mode_info{1}; mode_amps = mode_info{2};
    sumSize = length(tauAlphaN);

    a1_ = subplot(2,1,1);
    plot(1:sumSize, tauAlphaN, 'o-k', 'displayname',...
        'mode decay time constant')
    legend show
    set(a1_, 'xlim', [1 sumSize], 'xtick', 1:sumSize);
    xlabel('mode #')
    ylabel('mode time constant (ns)')
    grid

    a2_ = subplot(2,1,2);
    plot(tauAlphaN, mode_amps, 'o-k', 'displayname', ...
        'mode amplitude')
    xlabel('mode time constant (ns)')
    ylabel('mode amplitude (1/cm^3')
    hold all
    hD_ = findobj('tag','dataPlotLine'); 
    x = get(hD_, 'xdata'); y = get(hD_, 'ydata');
    plot(x, y, 'r.', 'displayname', 'data (incl shifts)')
    legend show
    grid
end

% Save the plotted model to the workspace:

plottedModel = struct;
plottedModel.x = t;
plottedModel.y = nAvg;
assignin('base', 'plottedModel', plottedModel);

end

% -------- Plots the data ----------------------------------------------
function plotData(handles)

axes(handles.mainAxes);

x = handles.dataX;
y = handles.dataY;

% Apply x shift, mobility scaling (if requested), and then the y shift.
% This will be plotted, recorded in the data table, and will be used when
% fitting. But the stored data will not be changed:

x = x + handles.dataXShift;
if get(handles.input_scaleByMobility, 'Value')
    y = y*handles.nomMajMob*handles.effMassRatio/handles.mu;
end
y = y + handles.dataYShift;

% Update the plot:

if length(x) == length(y)
    set(handles.text_dataError,'visible','off');
    hD_ = findobj('tag', 'dataPlotLine'); 
    if ~isempty(hD_)
        set(hD_, 'xdata', x,'ydata', y);
    else
        plot(x, y, 'r.', 'tag', 'dataPlotLine');
    end
else
    set(handles.text_dataError, 'visible', 'on');
end
grid on
box on

% Update the data table:

if ~iscolumn(x)
    x = x';
end
if ~iscolumn(y)
    y = y';
end
lX = length(x); lY = length(y);
if lX > lY
    xTable = x;
    yTable = [y; NaN(lX-lY,1)];
elseif lY > lX;
    xTable = [x; NaN(lY-lX,1)];
    yTable = y;
else
    xTable = x; yTable = y;
end
fDataTables_ = findobj('tag','fDataTables');
if isempty(fDataTables_)
    [tableModel_, tableData_] = makeTableFig;
else
    tableData_ = findobj('tag','tableData');
end
set(tableData_,'Data',[xTable yTable]);

% Write plotted data to workspace:

plottedData = struct;
plottedData.x = x;
plottedData.y = y;
assignin('base', 'plottedData', plottedData);

end


function [tableModel_, tableData_] = makeTableFig()

f_ = figure('Units','Pixels','Position',[70 90 472 370],...
    'tag','fDataTables','name','Number Tables');
tableModel_ = uitable('parent',f_,'tag','tableModel');
set(tableModel_,'ColumnName',{'t (ns)' 'n (um^-3)'},...
    'Units','Pixels','Position',[33 52 200 250],...
    'ColumnWidth',{89 89},...
    'RowStriping','on','RowName',[]);
tableData_ = uitable('parent',f_,'tag','tableData');
set(tableData_,'ColumnName',{'t (ns)' 'n (um^-3)'},...
    'Units','Pixels','Position',[250 52 200 250],...
    'ColumnWidth',{89 89},...
    'RowStriping','on','RowName',[]);
annotation(f_,'textbox',...
    [0.12 0.86 0.13 0.07],...
    'String',{'Model'});
annotation(f_,'textbox',...
    [0.67 0.87 0.11 0.07],...
    'String',{'Data (incl shifts)'},...
    'EdgeColor',[1 0 0]);

end


function input_SRV_Callback(hObject, eventdata, handles)
% hObject    handle to input_SRV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_SRV as text
%        str2double(get(hObject,'String')) returns contents of input_SRV as a double
handles.SRV = str2double(get(hObject,'String'));
updateInputs(handles);
plotModel(handles);
guidata(hObject,handles);
end


% --- Executes on slider movement.
function slider_SRV_Callback(hObject, eventdata, handles)
% hObject    handle to slider_SRV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.SRV = get(hObject,'Value');
updateInputs(handles);
plotModel(handles);
guidata(hObject,handles);
end

function input_difu_Callback(hObject, eventdata, handles)
% hObject    handle to input_difu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_difu as text
%        str2double(get(hObject,'String')) returns contents of input_difu as a double
handles.difu = str2double(get(hObject,'String'));
q = 1;
kB = 8.617e-5; % eV/K
handles.mu = handles.difu*q/kB/handles.temp;
updateInputs(handles);
plotModel(handles);
guidata(hObject,handles);
end

% --- Executes on slider movement.
function slider_difu_Callback(hObject, eventdata, handles)
% hObject    handle to slider_difu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.difu = get(hObject,'Value');
q = 1;
kB = 8.617e-5; % eV/K
handles.mu = handles.difu*q/kB/handles.temp;
updateInputs(handles);
plotModel(handles);
guidata(hObject,handles);
end


function input_tau_Callback(hObject, eventdata, handles)
% hObject    handle to input_tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_tau as text
%        str2double(get(hObject,'String')) returns contents of input_tau as a double
handles.tau = str2double(get(hObject, 'String'));
updateInputs(handles);
plotModel(handles);
guidata(hObject, handles);
end


% --- Executes on slider movement.
function slider_tau_Callback(hObject, eventdata, handles)
% hObject    handle to slider_tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.tau = get(hObject,'Value');
updateInputs(handles);
plotModel(handles);
guidata(hObject,handles);
end


function input_d_Callback(hObject, eventdata, handles)
% hObject    handle to input_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_d as text
%        str2double(get(hObject,'String')) returns contents of input_d as a double
handles.thick = str2double(get(hObject,'String'));
plotModel(handles);
guidata(hObject,handles);
end


function input_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to input_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_alpha as text
%        str2double(get(hObject,'String')) returns contents of input_alpha as a double
handles.alpha = str2double(get(hObject,'String'));
plotModel(handles);
guidata(hObject,handles);
end


function input_R_Callback(hObject, eventdata, handles)
% hObject    handle to input_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_R as text
%        str2double(get(hObject,'String')) returns contents of input_R as a double
handles.R = str2double(get(hObject,'String'));
plotModel(handles);
guidata(hObject,handles);
end


function input_mu_Callback(hObject, eventdata, handles)
% hObject    handle to input_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_mu as text
%        str2double(get(hObject,'String')) returns contents of input_mu as a double
q = 1;
kB = 8.617e-5; % eV/K
handles.mu = str2double(get(hObject,'String'));
handles.difu = handles.mu*kB*handles.temp/q;
set(handles.input_difu,'String',handles.difu);
updateInputs(handles);
plotModel(handles);
guidata(hObject,handles);
end


function input_temp_Callback(hObject, eventdata, handles)
% hObject    handle to input_temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_temp as text
%        str2double(get(hObject,'String')) returns contents of input_temp as a double
q = 1;
kB = 8.617e-5; % eV/K
handles.temp = str2double(get(hObject,'String'));
handles.difu = handles.mu*kB*handles.temp/q;
set(handles.input_difu,'String',handles.difu);
updateInputs(handles);
plotModel(handles);
guidata(hObject,handles);
end


% --- Executes during object creation, after setting all properties.
function menu_pumpProfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu_pumpProfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', {'Delta'; 'Square'; 'Gaussian'});
set(hObject, 'Value', 1)
end


% --- Executes on selection change in menu_pumpProfile.
function menu_pumpProfile_Callback(hObject, eventdata, handles)
% hObject    handle to menu_pumpProfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menu_pumpProfile contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menu_pumpProfile

menuStrings = get(hObject,'String');
handles.pumpProfile = menuStrings{get(hObject,'Value')};
plotModel(handles);
guidata(hObject,handles);

end


function input_N_Callback(hObject, eventdata, handles)
% hObject    handle to input_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_N as text
%        str2double(get(hObject,'String')) returns contents of input_N as a double

handles.N = str2double(get(hObject, 'String'));
plotModel(handles);
guidata(hObject, handles);

end


function input_T_Callback(hObject, eventdata, handles)
% hObject    handle to input_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_T as text
%        str2double(get(hObject,'String')) returns contents of input_T as a double

handles.T = str2double(get(hObject, 'String'));
plotModel(handles);
guidata(hObject, handles);

end


function input_h_Callback(hObject, eventdata, handles)
% hObject    handle to input_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_h as text
%        str2double(get(hObject,'String')) returns contents of input_h as a double

handles.h = str2double(get(hObject, 'String'));
plotModel(handles);
guidata(hObject, handles);

end


function input_sigma_Callback(hObject, eventdata, handles)
% hObject    handle to input_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_sigma as text
%        str2double(get(hObject,'String')) returns contents of input_sigma as a double

handles.sigma = str2double(get(hObject, 'String'));
plotModel(handles);
guidata(hObject, handles);

end


% --- Executes on button press in input_yAutoScale.
function input_yAutoScale_Callback(hObject, eventdata, handles)
% hObject    handle to input_yAutoScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of input_yAutoScale
as = get(hObject,'Value');
if as
    set(handles.mainAxes,'yLimMode','auto');
else
    set(handles.mainAxes,'yLimMode','manual');
end
updateInputs(handles);
end
        
function input_yAxesMin_Callback(hObject, eventdata, handles)
% hObject    handle to input_yAxesMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_yAxesMin as text
%        str2double(get(hObject,'String')) returns contents of input_yAxesMin as a double
yLim = get(handles.mainAxes,'YLim');
yLim(1) = str2double(get(hObject,'String'));
set(handles.mainAxes,'YLim',yLim);
set(handles.input_yAutoScale,'Value',0);
guidata(hObject, handles);
end

function input_yAxesMax_Callback(hObject, eventdata, handles)
% hObject    handle to input_yAxesMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_yAxesMax as text
%        str2double(get(hObject,'String')) returns contents of input_yAxesMax as a double
yLim = get(handles.mainAxes,'YLim');
yLim(2) = str2double(get(hObject,'String'));
set(handles.mainAxes,'YLim',yLim);
set(handles.input_yAutoScale,'Value',0);
guidata(hObject, handles);
end

% --- Executes on button press in input_xAutoScale.
function input_xAutoScale_Callback(hObject, eventdata, handles)
% hObject    handle to input_xAutoScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of input_xAutoScale
as = get(hObject,'Value');
if as
    set(handles.mainAxes,'xLimMode','auto');
else
    set(handles.mainAxes,'xLimMode','manual');
end
updateInputs(handles);
end

function input_xAxesMin_Callback(hObject, eventdata, handles)
% hObject    handle to input_xAxesMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_xAxesMin as text
%        str2double(get(hObject,'String')) returns contents of input_xAxesMin as a double
xLim = get(handles.mainAxes,'XLim');
xLim(1) = str2double(get(hObject,'String'));
set(handles.mainAxes,'XLim',xLim);
set(handles.input_xAutoScale,'Value',0);
guidata(hObject, handles);
end

function input_xAxesMax_Callback(hObject, eventdata, handles)
% hObject    handle to input_xAxesMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_xAxesMax as text
%        str2double(get(hObject,'String')) returns contents of input_xAxesMax as a double
xLim = get(handles.mainAxes,'XLim');
xLim(2) = str2double(get(hObject,'String'));
set(handles.mainAxes,'XLim',xLim);
set(handles.input_xAutoScale,'Value',0);
guidata(hObject, handles);
end


% --- Executes on button press in input_logScale.
function input_logScale_Callback(hObject, eventdata, handles)
% hObject    handle to input_logScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of input_logScale
sc = get(hObject,'Value');
if sc
    set(handles.mainAxes,'YScale','log');
else
    set(handles.mainAxes,'YScale','linear');
end
end


% --- Executes during object creation, after setting all properties.
function popupmenu_dataX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_dataX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','choose x')
end


% --- Executes on selection change in popupmenu_dataX.
function popupmenu_dataX_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_dataX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_dataX contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_dataX
vars = evalin('base','who');
ix = get(hObject,'Value')-1;
vX = vars{ix};
handles.dataX = evalin('base',vX);
plotData(handles)
guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function popupmenu_dataY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_dataY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','choose y')
end


% --- Executes on selection change in popupmenu_dataY.
function popupmenu_dataY_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_dataY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_dataY contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_dataY
vars = evalin('base','who');
iy = get(hObject,'Value')-1;
vY = vars{iy};
handles.dataY = evalin('base',vY);
plotData(handles)
guidata(hObject, handles);
end

function input_dataXShift_Callback(hObject, eventdata, handles)
% hObject    handle to input_dataXShift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_dataXShift as text
%        str2double(get(hObject,'String')) returns contents of input_dataXShift as a double
handles.dataXShift = str2double(get(hObject,'String'));
plotData(handles)
guidata(hObject, handles);
end

function input_dataYShift_Callback(hObject, eventdata, handles)
% hObject    handle to input_dataYShift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_dataYShift as text
%        str2double(get(hObject,'String')) returns contents of input_dataYShift as a double
handles.dataYShift = str2double(get(hObject,'String'));
plotData(handles)
guidata(hObject, handles);
end

% --- Executes on button press in input_scaleByMobility.
function input_scaleByMobility_Callback(hObject, eventdata, handles)
% hObject    handle to input_scaleByMobility (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of input_scaleByMobility
plotData(handles)
end


function input_effMassRatio_Callback(hObject, eventdata, handles)
% hObject    handle to input_effMassRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_effMassRatio as text
%        str2double(get(hObject,'String')) returns contents of input_effMassRatio as a double
handles.effMassRatio = str2double(get(hObject, 'String'));
plotData(handles)
plotModel(handles)
guidata(hObject, handles)
end


function input_nomMajMob_Callback(hObject, eventdata, handles)
% hObject    handle to input_nomMajMob (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_nomMajMob as text
%        str2double(get(hObject,'String')) returns contents of input_nomMajMob as a double
handles.nomMajMob = str2double(get(hObject, 'String'));
plotData(handles)
guidata(hObject, handles)
end



% --- Executes on button press in pushbutton_importData.
function pushbutton_importData_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_importData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vars = evalin('base','who');
set(handles.popupmenu_dataX,'String',['choose x';vars]);
set(handles.popupmenu_dataY,'String',['choose y';vars]);
set(handles.popupmenu_dataX,'Value',1);
set(handles.popupmenu_dataY,'Value',1);

handles.dataX = [];
handles.dataY = [];
guidata(hObject,handles);

delete(findobj('tag','dataPlotLine'));

tableData_ = findobj('tag','tableData');
if ~isempty(tableData_)
    set(tableData_,'Data',[]);
end

end


% --- Executes on button press in pushbutton_fit.
function pushbutton_fit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PNames = {'tau' 'SRV' 'thickness' 'alpha' 'R' 'diffusivity' 'N' 'sigma' 'T' ...
    'h' 'dataXShift' 'dataYShift'};
% array of names

q = 1;
kB = 8.617e-5; % eV/K
if get(handles.input_scaleByMobility, 'Value')
    dataScale = {true;...
        handles.nomMajMob*handles.effMassRatio*kB*handles.temp/q};
else
    dataScale = {false;...
        handles.nomMajMob*handles.effMassRatio*kB*handles.temp/q};
end
% data partial scaling factor. Needs to be divided by the minority
% carrier diffusivity, which is done inside the kernel in fitFC

P = [handles.tau, handles.SRV, handles.thick, handles.alpha, handles.R,...
    handles.difu, handles.N, handles.sigma, handles.T, handles.h, handles.dataXShift, ...
    handles.dataYShift];
% array of model parameters

fit = [get(handles.input_fit_tau,'Value'), get(handles.input_fit_SRV,'Value'),...
    0, get(handles.input_fit_alpha,'Value'), 0, get(handles.input_fit_difu,'Value'), ...
    get(handles.input_fit_N,'Value'), 0, 0, get(handles.input_fit_h, 'Value'),...
    get(handles.input_fit_dataXShift,'Value'), get(handles.input_fit_dataYShift,'Value')];
% array of logicals indicating which parameters to fit
if strcmp(handles.pumpProfile, 'Gaussian')
    fit(8) = get(handles.input_fit_sigma, 'Value');
end
fix = ~fit;
% logical arrays of what to fit, and what to fix       

xLim = get(handles.mainAxes, 'XLim');
x = handles.dataX; y = handles.dataY;
g = (x + handles.dataXShift >= min(xLim)) & (x + handles.dataXShift <= max(xLim));
% filter data based on current plot, including the X shift.

if get(handles.input_weights,'Value')
    xStep = zeros(size(x(g)));
    xStep(2:end) = diff(x(g));
    xStep(1) = xStep(2);
    w = xStep/min(xStep);
    % use normalized stepsizes as weights
    w = w.^2;
else
    w = ones(size(x(g)));
end


injDep = struct;
% structure to pass injection dependence info to solveFC
injDep.tauModel     = handles.tauModel;
injDep.difuModel    = handles.difuModel;
injDep.effMassRatio = handles.effMassRatio;
injDep.majConc      = handles.majConc;
injDep.temp         = handles.temp;

wait_ = [];
f_ = get(handles.mainAxes, 'parent');
wait_ = annotation(f_, 'textbox', [.15 .5 .2 .1], 'String',...
    'Fitting ...', 'EdgeColor', 'r', 'Color', 'r', 'BackgroundColor', 'w',...
    'FontSize', 25);
pause(.05)
% throw up a "wait" signal

fitOut = fitFC(x(g), y(g), w, handles.pumpProfile, P, fix, dataScale,...
    handles.solver, injDep);

if ~isempty(wait_)
    delete(wait_)
end

mdl     = fitOut.mdl;
POut    = fitOut.POut;
SEOut   = fitOut.SEOut;
PFit    = fitOut.PFit;
Fitted  = fitOut.Fitted;

% Update all the params in the GUI:
handles.tau         = POut(1); 
handles.SRV         = POut(2); 
handles.thick       = POut(3);
handles.alpha       = POut(4); 
handles.R           = POut(5); 
handles.difu        = POut(6);
handles.N           = POut(7); 
handles.sigma       = POut(8); 
handles.T           = POut(9);
handles.h           = POut(10); 
handles.dataXShift  = POut(11); 
handles.dataYShift  = POut(12);
q = 1; kB = 8.617e-5; % eV/K
handles.mu          = handles.difu*q/kB/handles.temp;

updateInputs(handles);
plotModel(handles);
plotData(handles);
guidata(hObject,handles);

% Update the fit number, save fit data to the base workspace:
fitNum = str2num(get(handles.indicate_fitNum, 'String'));
fitNum = fitNum + 1;
set(handles.indicate_fitNum, 'String', num2str(fitNum, '%u'));

fitInfo = struct;
% structure to save information about the fit; used subsequently for
% various things
fitInfo.model           = mdl;
% the fit model returned by fitnlm() !!!!!!!!!!!!!
fitInfo.datasetName     = get(handles.input_datasetName, 'String');
fitInfo.x               = x(g);
fitInfo.y               = y(g);
fitInfo.w               = w;
% the data that was fit, and the weights
fitInfo.Fitted          = Fitted;
% the model evaluated at the best-fit parameters
fitInfo.pumpProfile     = handles.pumpProfile;
% the pump profile
fitInfo.PFit            = PFit;
% initial values of the fitted parameters
fitInfo.POut            = POut;
% best-fit parameters
fitInfo.SEOut           = SEOut;
% standard error
fitInfo.PNames          = PNames;
% parameter names
fitInfo.PFitNames       = PNames(~fix);
% names of fitted parameters
fitInfo.fix             = fix;
% boolean array indicating the parameters that were held fixed
fitInfo.dataScale       = dataScale;
% structure indicating whether data is scaled by mobilit, and the nominal
% scaling factor
fitInfo.solver          = handles.solver;
% indicates the solver used ('Analytical' or 'Numerical')
fitInfo.injDep          = injDep;
% structure describing the model to describe injection-dependence of the
% recombination rates

assignin('base', sprintf('fitInfo%u', fitNum), fitInfo);
% save fit info structure to the base workspace

% Report results in a table:

fLastModel_ = findobj('tag', 'fLastModel'); delete(fLastModel_);
fLastModel_ = figure('name', 'Last Fit Model', 'tag', 'fLastModel',...
    'Position', [900 45 415 360]);
tModelParams_ = uitable('parent', fLastModel_, 'tag', 'tModelParams');
set(tModelParams_, 'ColumnName', {'Parameter' 'Fit?' 'Value' '95% CI'});
set(tModelParams_, 'ColumnFormat', {'char' 'logical' 'numeric' 'numeric'});
set(tModelParams_, 'Units', 'Pixels', 'Position', [45 20 354 266]);
set(tModelParams_, 'FontSize', 10, 'RowName', []);
set(tModelParams_, 'ColumnWidth', {100 50 100 100});

rowNames = {'tau (ns)'; 'SRV (cm/s)'; 'thick (um)'; 'alpha (1/cm)'; 'R (0-1)';...
    'diffusivity (cm^2/s)'; 'N (1/cm^2)'; 'sigma (ns)'; 'T (ns)'; 'h (1/cm)'; ...
    'xShift (ns)'; 'yShift (1/um^3)'};
colData = cell(length(rowNames), 1);
for i = 1:length(rowNames)
   colData{i,1} = rowNames{i};
   colData{i,2} = logical(fit(i));
   colData{i,3} = POut(i);
   colData{i,4} = SEOut(i);
   if ~fit(i) & isnan(SEOut(i))
       colData{i,4} = [];
   end
   % If the param was not fit, just leave CI entry blank rather than NaN
end
switch handles.pumpProfile
    case 'Delta'
        colData{8, 3} = [];
        colData{8, 4} = [];
        colData{9, 3} = [];
        colData{9, 4} = [];
    case 'Square'
        colData{8, 3} = [];
        colData{8, 4} = [];
end
% for delta and square pulses, don't show anything for the unused
% parameters
set(tModelParams_, 'Data', colData);
ann_ = annotation(fLastModel_, 'textbox','Position', [.11 .81 .87 .17], ...
    'EdgeColor', 'none', 'FontSize', 12, 'tag', 'pump annotation',...
    'Units', 'Normalized', 'String',...
    ['pump = ' handles.pumpProfile ', solver = ' handles.solver char(10)...
    'scale data by mobility = ' num2str(dataScale{1}) char(10)...
    'tauModel = ' handles.tauModel ', difuModel = ' handles.difuModel]);

% Plots to visualize covariance matrix:

if get(handles.input_plotCovariance, 'Value')
    covarPairPlots(fitInfo);
end

end


% --- Executes during object creation, after setting all properties.
function menu_tauModel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu_tauModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', {'Constant'; 'Simple SRH'});
set(hObject, 'Value', 1);
end


% --- Executes on selection change in menu_tauModel.
function menu_tauModel_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tauModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menu_tauModel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menu_tauModel

tauModelOptions = get(hObject, 'String');
if strcmp(handles.solver, 'Analytical')
    handles.tauModel = tauModelOptions{1};
    set(hObject, 'Value', 1)
else
    handles.tauModel = tauModelOptions{get(hObject, 'Value')};
end
plotModel(handles)
guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function menu_difuModel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu_difuModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', {'Constant'; 'Ambipolar'});
set(hObject, 'Value', 1);
end


% --- Executes on selection change in menu_difuModel.
function menu_difuModel_Callback(hObject, eventdata, handles)
% hObject    handle to menu_difuModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menu_difuModel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menu_difuModel

difuModelOptions = get(hObject, 'String');
if strcmp(handles.solver, 'Analytical')
    handles.difuModel = difuModelOptions{1};
    set(hObject, 'Value', 1)
else
    handles.difuModel = difuModelOptions{get(hObject, 'Value')};
end
handles.difuModel = difuModelOptions{get(hObject, 'Value')};
plotModel(handles)
guidata(hObject, handles);
end


function input_majConc_Callback(hObject, eventdata, handles)
% hObject    handle to input_majConc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_majConc as text
%        str2double(get(hObject,'String')) returns contents of input_majConc as a double

handles.majConc = str2double(get(hObject, 'String'));
plotModel(handles)
guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function menu_solver_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu_solver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', {'Analytical'; 'Numerical'});
set(hObject, 'Value', 1);
end


% --- Executes on selection change in menu_solver.
function menu_solver_Callback(hObject, eventdata, handles)
% hObject    handle to menu_solver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menu_solver contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menu_solver

solverOptions = get(hObject, 'String');
handles.solver = solverOptions{get(hObject, 'Value')};
if strcmp(handles.solver, 'Analytical')
   handles.tauModel = 'Constant';
   handles.difuModel = 'Constant';
end
updateInputs(handles)
plotModel(handles)
guidata(hObject, handles);

end


% --- Executes on button press in input_latchModel.
function input_latchModel_Callback(hObject, eventdata, handles)
% hObject    handle to input_latchModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of input_latchModel
if ~get(hObject, 'Value')
    plotModel(handles);
end
end
