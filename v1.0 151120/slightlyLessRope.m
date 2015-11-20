function slightlyLessRope(fitInfoArray)

% function slightlyLessRope(fitInfoArray)
%
% Simultaneously fits multiple datasets from transiet
% free-carrier-absorption experiments. The data manipulation and fit models
% are as in plentyOfRope.m. The input is an array of fit information
% structures as produced by plentyOfRope.m. The can then select how to
% fit and constrain varialbles and then perform a global fit accross all of
% the datasets.
%
% fitInfoArray = array of fit information structures as saved by
% plentyOfRope.m.
%
% July 2015, R Jaramillo: modified from flexFitMultiLukeAvgFC.m.

N = length(fitInfoArray);

solverChoices       = {'Analytical' 'Numerical'};
tauModelChoices     = {'Constant' 'Simple SRH'};
difuModelChoices    = {'Constant' 'Ambipolar'};

paramLabels = {  'tau (ns)'; 
                'SRV (cm/s)'; 
                'thickness (um)'; 
                'alpha (1/cm)'; 
                'R (0-1)';
                'D (cm^2/s)'; 
                'N (1/cm^2)'; 
                'sigma (ns)'; 
                'T (ns)'; 
                'h (1/cm)'; 
                'xShift (ns)'; 
                'yShift (1/um^3)'};
colNames = {    'value';
                'was fit?';
                'fit now?';
                'fit value'}; 
colFormat = {   'numeric'...
                'logical'...
                'logical'};
paramTags = {'tau';
                'S';
                'th';
                'a';
                'R';
                'D';
                'N';
                'sig';
                'T';
                'h';
                'xS';
                'yS';};
colEditable = [0 0 1 1];

allColNames     = {'dataset name'; 'Rescaled by |mobility?'; 'Rescale |now?';...
    'Solver used?'; 'Solver now?'; 'tauModel used?'; 'tauModel now?';...
    'difuModel used?'; 'difuModel now?'};
allColFormat    = {'char' 'logical' 'logical' solverChoices solverChoices ...
    tauModelChoices tauModelChoices difuModelChoices difuModelChoices};
allColEditable  = [0 0 1 0 1 0 1 0 1];
allColWidth     = {350 'auto' 'auto' 'auto' 90 'auto' 110 'auto' 'auto'};
allDat          = {};

for j = 1:length(paramLabels)
    
    l = paramTags{j};
    varNames = {};
    for k = 1:20
       varNames = [varNames, [l sprintf('%u', k)]]; 
    end
    
    colNames_j = colNames;
    for k = 1:length(colNames)
       colNames_j{k} = [paramLabels{j} ' |' colNames_j{k}]; 
    end
    
    allColNames     = [allColNames;     colNames_j];
    allColFormat    = [allColFormat,    [colFormat {varNames}]];
    allColEditable  = [allColEditable,  colEditable];
    allColWidth     = [allColWidth, 85, 'auto', 'auto', 'auto'];
    
end

for j = 1:N
    
    POut            = fitInfoArray(j).POut;
    fix             = fitInfoArray(j).fix;
    dataScale       = fitInfoArray(j).dataScale;    
    nextDatasetName = [fitInfoArray(j).datasetName ', ' fitInfoArray(j).pumpProfile];   
    solver          = fitInfoArray(j).solver;
    injDep          = fitInfoArray(j).injDep;
    tauModel        = injDep.tauModel;
    difuModel       = injDep.difuModel;

    nextDat = cell( 1, 9 + 4*length(paramLabels) );
    
    nextDat{1}  = nextDatasetName;
    nextDat{2}  = dataScale{1};
    nextDat{3}  = dataScale{1};
    nextDat{4}  = solver;
    nextDat{5}  = solver;
    nextDat{6}  = tauModel;
    nextDat{7}  = tauModel;
    nextDat{8}  = difuModel;
    nextDat{9}  = difuModel;
   
   for k = 1:length(paramLabels)
   
      nextDat{ 1, 9 + 4*(k-1) + 1 } = POut(k);
      nextDat{ 1, 9 + 4*(k-1) + 2 } = ~fix(k);
      nextDat{ 1, 9 + 4*(k-1) + 3 } = false;
      nextDat{ 1, 9 + 4*(k-1) + 4 } = '';
   
   end
   
   allDat = [allDat; nextDat];
   
end


f_ = figure('Position', [610 160 1350 460], 'name', 'slightlyLessRope');
guidata(f_, fitInfoArray);
% store the main function input as the gui data

tDat_   = uitable(f_, 'ColumnName', allColNames, 'ColumnFormat', allColFormat,...
    'ColumnEditable', logical(allColEditable), 'RowName', 'numbered', ...
    'CellEditCallback', @table_CellEditCallback, 'tag', 'dataTable',...
    'ColumnWidth', allColWidth, 'Position', [20 20 1300 330],...
    'FontSize', 11);
set(tDat_, 'Data', allDat);

fit_ = uicontrol(f_, 'Style', 'pushbutton', 'String', 'Do Fit',...
    'Callback', @doFlexFit_Callback, 'Position', [20 390 80 60],...
    'FontSize', 12, 'BackgroundColor', [0.04 0.5 0.8]);

plotCov_ = uicontrol(f_, 'Style', 'checkbox', 'String', 'Plot Covar?',...
    'Position', [130 390 20 20], 'Value', 0, 'tag', 'plotCov?');

plotCov_label_ = uicontrol(f_, 'Style', 'text', 'String', 'Plot Covar?', ...
    'Position', [145 386 75 20], 'HorizontalAlignment', 'left');

scaleWeights_ = uicontrol(f_, 'Style', 'checkbox', 'String', 'Scale Weights?',...
    'Position', [330 390 20 20], 'Value', 0, 'tag', 'scaleWeights?');

scaleWeights_label_ = uicontrol(f_, 'Style', 'text', 'String', 'Scale Weights?', ...
    'Position', [354 387 90 20], 'HorizontalAlignment', 'left');

end

%%
function doFlexFit_Callback(hObject, eventdata)

fitInfoArray    = guidata(hObject);
N               = length(fitInfoArray);
tDat_           = findobj('tag', 'dataTable');
tDat            = get(tDat_, 'Data');
PNames          = fitInfoArray(1).PNames;
% parameter names for a given fit (should be the same for all fits)

lb = [0 0 0 0 0 0 0 0 -Inf 0 -Inf -Inf];
ub = [Inf Inf Inf Inf 1 Inf Inf Inf Inf Inf Inf Inf];
% lower and upper bounds on parameters 

xIn             = [];
yIn             = [];
wIn             = [];
PFitTags        = {};
% letter indices of fit variables 
PFitValues      = [];
% starting (and final) values of fit variables
PFitLB          = [];
PFitUB          = [];
% upper and lower bounds on fit variables
PFitNames       = {};
% names of fitted parameters

scW = get(findobj('tag', 'scaleWeights?'), 'Value');


for i = 1:N
    
    if ~all(strcmp(fitInfoArray(i).PNames, PNames))
        disp('dummy dummy dummy')
    end
    % make sure that every function has the same parameter names
    
    xNow = fitInfoArray(i).x;
    yNow = fitInfoArray(i).y;
    wNow = fitInfoArray(i).w;
    
    if scW
       maxY = max(yNow); 
       wNow = wNow/maxY^2;
    end
    % scale weights by the maximum value of the data. Ensures that scans
    % with different scales are counted equally in the chi-squared routine

    if ~isrow(xNow)
        xNow = xNow';
    end
    if ~isrow(yNow)
        yNow = yNow';
    end
    if ~isrow(wNow)
        wNow = wNow';
    end
    % make sure data are row vectors
    
    xIn = [xIn xNow];
    yIn = [yIn yNow];
    wIn = [wIn wNow];
    % append all of the x, y, and w datasets into long arrays
    
    for j = 1:length(PNames)
        if (sum(strcmp(tDat{i, 9 + 4*j}, PFitTags)) == 0) ...
                & (~isempty(tDat{i, 9 + 4*j}))
            PFitTags    = [PFitTags tDat{i, 9 + 4*j}];
            PFitValues  = [PFitValues tDat{i, 9 + 4*j - 3}];
            PFitLB      = [PFitLB lb(j)];
            PFitUB      = [PFitUB ub(j)];
            name_j = [PNames{j} ' (variable ''' tDat{i, 9 + 4*j} ''')'];
            PFitNames = [PFitNames name_j];
        end
    end
    % populate arrays representing the variables to be fit
    
end

if length(PFitValues)==0
    disp('Nothing to fit you silly rabbit')
    return
end

PFitNorm    = ones(size(PFitValues));

% opts        = statset('fitnlm'); opts = statset(opts, 'Display', 'iter');
% mdl         = fitnlm(xIn, yIn, @fun, PFitNorm, 'Weights', wIn, 'Options', opts);
% PFitOut     = (mdl.Coefficients.Estimate').*PFitValues;
% SEOut       = (mdl.Coefficients.SE').*abs(PFitValues);

opts = optimoptions('lsqcurvefit', 'Display', 'iter-detailed');
mdl = struct;
[mdl.x, mdl.resnorm, mdl.residual, mdl.exitflag, mdl.output, mdl.lambda,...
    mdl.jacobian] = lsqcurvefit(@fun, PFitNorm, xIn, yIn.*sqrt(wIn), PFitLB, PFitUB, opts);
ci          = nlparci(mdl.x, mdl.residual, 'jacobian', mdl.jacobian);
PFitOut     = (mdl.x).*PFitValues;
SEOut       = abs(ci(:, 1)' - mdl.x).*abs(PFitValues);

% Create output structure:
fitOut                  = struct;
fitOut.fitInfoArray     = fitInfoArray;
fitOut.model            = mdl;
fitOut.PFit             = PFitValues;
fitOut.PFitOut          = PFitOut;
fitOut.SEOut            = SEOut;
fitOut.Fitted           = (mdl.residual)./sqrt(wIn) + yIn;
fitOut.PFitNames        = PFitNames;

%% Create figure and table to report results:
delete(findobj('tag', 'multiFitResult'));
delete(findobj('tag', 'covarResult'));
% delete output figures from previous fits

fOut_ = figure('Position', [150 150 1092 503], 'Name', 'slightlyLessRope_out',...
    'tag', 'multiFitResult', 'UserData', fitOut);
aOut_ = axes('Position', [0.047 0.11 0.269 0.815]);
hold(aOut_,'all');

paramLabels = {  'tau (ns)'; 
                'SRV (cm/s)'; 
                'thickness (um)'; 
                'alpha (1/cm)'; 
                'R (0-1)';
                'D (cm^2/s)'; 
                'N (1/cm^2)'; 
                'sigma (ns)'; 
                'T (ns)'; 
                'h (1/cm)'; 
                'xShift (ns)'; 
                'yShift (1/um^3)'};
colNames = {    'was fit?';
                'fit value';
                'value';
                '95% CI'}; 
colFormat = {   'logical'...
                'char'...
                'numeric'...
                'numeric'};

datOut = {};

allColNames     = {'dataset name'; 'Recaled by mobility?'; 'Solver used?';...
    'tauModel used?'; 'difuModel used?'};
allColFormat    = {'char' 'logical' 'char' 'char' 'char'};
allColWidth     = {350, 'auto' 'auto' 'auto' 'auto'};

for j = 1:length(paramLabels)
    
    colNames_j = colNames;
    for k = 1:length(colNames)
       colNames_j{k} = [paramLabels{j} ' |' colNames_j{k}]; 
    end
    
    allColNames     = [allColNames;     colNames_j];
    allColFormat    = [allColFormat,    colFormat];
    allColWidth     = [allColWidth, 'auto', 'auto', 85, 85];
    
end

startInd = 1;
% starting index of this dataset in the long arrays passed to and from
% lsqcurvefit()
yModel = fitOut.Fitted;
for i = 1:N
    
    x           = fitInfoArray(i).x;
    y           = fitInfoArray(i).y;
    P           = zeros(size(PNames));
    SE          = NaN(size(PNames));
    nextDat     = cell(1, 5 + 4*length(PNames));
    nextDat{1}  = tDat{i, 1};
    nextDat{2}  = tDat{i, 3};
    nextDat{3}  = tDat{i, 5};
    nextDat{4}  = tDat{i, 7};
    nextDat{5}  = tDat{i, 9};
    
    stopInd = startInd + length(x) - 1;
    % starting index of this dataset in the long arrays passed to and from
    % fitnlm()
    
    for j = 1:length(PNames)
        
        P(j)    = tDat{i, 9 + 4*j - 3};
        g = strcmp(tDat{i, 9 + 4*j}, PFitTags);
        if sum(g) == 1
            P(j)    = PFitOut(g);
            SE(j)   = SEOut(g);
        elseif sum(g) > 1
            disp('dummy dummy dummy dummy dummy');
            return
        end
        
        nextDat{5 + 4*j - 3} = tDat{i, 9 + 4*j - 1};
        nextDat{5 + 4*j - 2} = tDat{i, 9 + 4*j};
        nextDat{5 + 4*j - 1} = P(j);
        nextDat{5 + 4*j}     = SE(j);
        
    end
    
    datOut = [datOut; nextDat];
    % populate arrays of parameters P and 95% confidence interval ("SE"), and build the
    % data for the output table
    
    
    yModel_i = yModel(startInd:stopInd);
    
    dataXShift = P(end-1);
    dataYShift = P(end);
    
    x = x + dataXShift;
    if tDat{i, 3}
        dataScale   = fitInfoArray(i).dataScale;
        scFun       = dataScale{2}/P(6);
        y           = y*scFun + dataYShift;
        yModel_i    = yModel_i*scFun + dataYShift;
    else
        y           = y + dataYShift;
        yModel_i    = yModel_i + dataYShift;
    end
    l_ = plot(x, y, '.', 'parent', aOut_,...
        'displayname', fitInfoArray(i).datasetName);
    plot(x, yModel_i, '-', 'displayname', 'model', 'color', ...
        get(l_, 'color'));
    
    startInd = stopInd + 1;
end

xlabel('time (ns)')
ylabel('\Deltan (cm^{-3})')
box on
legend show

tOut_   = uitable(fOut_, 'ColumnName', allColNames, 'ColumnFormat', allColFormat,...
    'ColumnEditable', [], 'RowName', 'numbered', 'Position', [370 18 716 321],...
    'FontSize', 11, 'ColumnWidth', allColWidth);
set(tOut_, 'Data', datOut);

if scW
    scWAnn = 'ON';
else
    scWAnn = 'OFF';
end
ann_ = annotation(fOut_, 'textbox', 'Units', 'Normalized', 'Position', ...
    [.5 .75 .49 .15], 'EdgeColor', [0 .5 0], 'FontSize', 12, 'tag', ...
    'BestFitAnnotation', 'Color', 'r', 'String', ...
    ['scale Weights? = ' scWAnn]);

% Create covariance pair plots
if get(findobj('tag', 'plotCov?'), 'Value')
    fCov_ = covarPairPlots(fitOut);
    set(fCov_, 'tag', 'covarResult');
end

%% Function for fitting routine
    function f = fun(PFun, xFun)
        % fun() = function to be fit by fitnlm()
        % PFun = parameters being fit. All normalized by the starting
        % values.
        % xFun = passed but not used.
        
        fTmp = [];
        fTmp_noScale = [];
        
        for j = 1:N
            
            P = [];
            for k = 1:length(PNames)
               P = [P; tDat{j, 9 + 4*k - 3}];
            end
            % starting values of all parameters (not normalized)
            
            for k = 1:length(P)
                g = strcmp(tDat{j, 9 + 4*k}, PFitTags);
                if sum(g) == 1
                    P(k) = PFun(g)*PFitValues(g);
                elseif sum(g) > 1
                    disp('dummy dummy dummy dummy')
                    return
                end
            end
            % step through the parameters and update the values being fit
            
            PSolveFC = P(1:9);
            % parameters that are passed to solveFC (all but xShift and
            % yShift)
            h = P(10);
            xShift = P(11);
            yShift = P(12);      
            pump                = fitInfoArray(j).pumpProfile;
            solver              = tDat{j, 5};
            injDep              = fitInfoArray(j).injDep;
            injDep.tauModel     = tDat{j, 7};
            injDep.difuModel    = tDat{j, 9};
            % update the solver and injDep structure to reflect choices in GUI
            
            xNow = fitInfoArray(j).x;
            yNow = fitInfoArray(j).y;
            wNow = fitInfoArray(j).w;
            if scW
               maxY = max(yNow); 
               wNow = wNow/maxY^2;
            end
            if ~isrow(xNow)
                xNow = xNow';
            end
            if ~isrow(wNow)
                wNow = wNow';
            end
            [y, ~, ~] = solveFC(xNow + xShift, NaN, NaN, pump, PSolveFC, h,...
                solver, injDep, false);
            
            dataScale = fitInfoArray(j).dataScale;
            if tDat{j, 3}
                scFun = dataScale{2}/P(6);
                % scale function, depends on diffusivity
                y = (y - yShift)/scFun;;
                % apply the y-shift and scaling
            else
                y = y - yShift;
            end
            
            fTmp_noScale = [fTmp_noScale y];
            
            y = y.*sqrt(wNow);
            % apply weights (also applied inside lsqnonline.m function call)
            
            fTmp = [fTmp y];

        end
        
        f = fTmp;
        if (isrow(xFun) & ~isrow(f)) | (iscolumn(xFun) & ~iscolumn(f))
            f = f';
        end
        % makes sure that input and output are the same size
    
    end

end


% --- Executes when entered data in editable cell(s) in the data table
function table_CellEditCallback(hObject, eventdata)
% hObject    handle to the table that was edited (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data

% This callback ensures that if a parameter is not being fit, then nothing
% can be entered in the 'fit value' column. Also ensures that only relevant
% choices can be entered for 'tauModel now?' and 'diffModel now?'.

ind = eventdata.Indices;
dat = get(hObject, 'Data');

if ind(2) > 9 & mod(ind(2)-9, 4) == 0
    if ~dat{ind(1), ind(2)-1}
        dat{ind(1), ind(2)} = '';
        set(hObject, 'Data', dat);
    end
end
% If user entered a fit value, and 'fit now?' is false, then remove the
% entry

if ind(2) > 9 & mod(ind(2)-9, 3) == 0
    if ~dat{ind(1), ind(2)}
        dat{ind(1), ind(2)+1} = '';
        set(hObject, 'Data', dat);
    end
end
% If user changes 'fit now?' to false, then remove the fit value entry

if ind(2) == 5
   if strcmp(dat{ind(1), ind(2)}, 'Analytical')
       dat{ind(1), 7} = 'Constant';
       dat{ind(1), 9} = 'Constant';
       set(hObject, 'Data', dat);
   end
    
end
% If user changed solver to Analytical, make sure that tauModel and
% difuModel are both Constant

if ind(2) == 7
    if strcmp(dat{ind(1), 5}, 'Analytical')
        dat{ind(1), ind(2)} = 'Constant';
        set(hObject, 'Data', dat);
    end 
end
% Make sure that 'tauModel now?' is 'Constant' in case the 'Solver now?'
% field is 'Analytical'

if ind(2) == 9
    if strcmp(dat{ind(1), 5}, 'Analytical')
        dat{ind(1), ind(2)} = 'Constant';
        set(hObject, 'Data', dat);
    end 
end
% Make sure that 'difuModel now?' is 'Constant' in case the 'Solver now?'
% field is 'Analytical'

end

