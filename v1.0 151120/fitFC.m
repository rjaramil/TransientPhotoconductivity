function fitOut = fitFC(xIn, yIn, wIn, pump, P, fix, dataScale, solver, injDep)

% function fitOut = fitFC(xIn, yIn, wIn, pump, P, fix, dataScale, solver, injDep)
%
% Fits a model of transient free carrier dynamics to data. Can use either the
% analytic model of Luke & Cheng (1987), or a numerical solver that allows
% for nonlinear terms in the model (e.g. high injection)
%
% xIn, yIn = the data to be fit. Not including X and Y shifts, and Y
% scaling.
%
% wIn = an array of weights to use during fitting. Should be same
% size as xIn and yIn
%
% pump = specifies optical pump profile; options are 'Delta', 'Square', and 'Gaussian' 
%
% P = model parameters. Exact parameter set depends on pump
%       Delta:      P = [tau, SRV, thick, alpha, R, difu, N, ~,     ~,  h, xShift, yShift]
%       Square:     P = [tau, SRV, thick, alpha, R, difu, N, ~,     T,  h, xShift, yShift]
%       Gaussian:   P = [tau, SRV, thick, alpha, R, difu, N, sigma, T,  h, xShift, yShift]
%
% fix = array of logicals (same size as PIn) specifying which parameters to
% hold fixed
%
% dataScale = indicates whether the data should be scaled by the mobility. Should
% be a cell array.
%   dataScale{1} = logical, indicating whether data should be scaled
%   dataScale{2} = the partial scale factor, to be divided by the minority carrier
%   diffusivity 
%       = (nominal majority carrier mobility) 
%       x (ratio of majority/minority carrier effective masses)
%
% solver = selects the solver: 'Analytical', or 'Numerical'. Passed
%   directly to solveFC.m
%
% injDep = structure to describe the injection-dependence of the model.
% Only applies to Numerical solver. Passed directly to solveFC.m
%   injDep.tauModel = selects bulk lifetime model
%   injDep.difuModel = selects diffusion model
%   injDep.effMassRatio = majority/minority carrier effective mass ratio
%   injDep.majConc = majority carrier concentration (cm^-3)
%
% fitOut = structure containing:
%   POut = best-fit parameters
%   SEOut = 95% confidence interval on fitted parameters
%   mdl = nonlinear model 
%
% July 2015, R. Jaramillo

lb = [0 0 0 0 0 0 0 0 -Inf 0 -Inf -Inf];
ub = [Inf Inf Inf Inf 1 Inf Inf Inf Inf Inf Inf Inf];
% lower and upper bounds on parameters

if ~isrow(xIn)
    xIn = xIn';
end
if ~isrow(yIn)
    yIn = yIn';
end
if ~isrow(wIn)
    wIn = wIn';
end
% make sure input data are row vectors


PFit = P(~fix);
PFitNorm = ones(size(PFit)); % Scale fit parameters; start with all values = 1

% opts = statset('fitnlm'); opts = statset(opts, 'Display', 'iter');
% mdl = fitnlm(xIn, yIn, @fun, PFitNorm, 'Weights', wIn, 'Options', opts);

% POut = P;
% POut(~fix) = (mdl.Coefficients.Estimate').*PFit;
% SEOut = NaN(size(POut));
% SEOut(~fix) = (mdl.Coefficients.SE').*abs(PFit);


opts = optimoptions('lsqcurvefit', 'Display', 'iter-detailed');
mdl = struct;
[mdl.x, mdl.resnorm, mdl.residual, mdl.exitflag, mdl.output, mdl.lambda,...
    mdl.jacobian] = lsqcurvefit(@fun, PFitNorm, xIn, yIn.*sqrt(wIn), lb(~fix), ub(~fix), opts);
ci = nlparci(mdl.x, mdl.residual, 'jacobian', mdl.jacobian);
POut = P;
POut(~fix) = (mdl.x).*PFit;
SEOut = NaN(size(POut));
SEOut(~fix) = abs(ci(:, 1)' - mdl.x).*abs(PFit);

% Create the output structure:
fitOut          = struct;
fitOut.POut     = POut;
fitOut.SEOut    = SEOut;
fitOut.mdl      = mdl;
fitOut.PFit     = PFit;
fitOut.Fitted   = mdl.residual + yIn;

% The model function:
    function f = fun(PFun, xFun)
        % fun() = function to be fit by fitnlm()
        % PFun = parameters being fit. In general will be a subset of the
        % full parameter set P. Normalized by the starting values.
        
        PFull = P;
        % full parameter set
        
        PFull(~fix) = PFun.*PFit;
        % update the parameters being fit. Rescale parameters.
        
        PSolveFC = PFull(1:9);
        % model parameters calculated by LukeAvgFC() or numFC()
        
        h = PFull(10);
        xShift = PFull(11);
        yShift = PFull(12);
        
        % Evaluate the model, with the x-data shifted:
        [f, ~, ~] = solveFC(xFun + xShift, NaN, NaN, pump,...
            PSolveFC, h, solver, injDep, false);
                
        if dataScale{1}
            scFun = dataScale{2}/PFull(6);
            % scale function, depends on diffusivity
            f = (f - yShift)/scFun;
            % apply the y-shift and scaling
        else
            f = f - yShift;
            % apply just the y-shift
        end
        
        f = f.*sqrt(wIn);
        
        if (isrow(xFun) & ~isrow(f)) | (iscolumn(xFun) & ~iscolumn(f))
            f = f';
        end
        % makes sure that input and output are the same size
        
    end

end

