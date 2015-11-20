function f_ = covarPairPlots(fitInfo)

% function f_ = covarPairPlots(fitInfo)
%
% Creates graphical representations of the covariance matrix for a
% nonlinear regression. Assumes that all fit parameters were normalized to
% unity before the fit.
%
% fitInfo = fit info structure saved by plentyOfRope.m or
% slightlyLessRope.m
%
% RJ, fall 2014
%
% Adapted to take input structure from plentyOfRope.m. July 2015, R
% Jaramillo
%
% Adapted to take input from lsqcurvefit() instead of fitnlm(), July 2015,
% R. Jaramillo

mdl         = fitInfo.model;
% the fit model
PFit        = fitInfo.PFit;
% initial values of the fitted parameters
PFitNames   = fitInfo.PFitNames;
% cell array of parameter names

N = length(mdl.x);
if N < 2
    f_ = [];
    return
end
M = nchoosek(N,2); % number of pairwise comparisons; number of plots
L = ceil(sqrt(M)); % the smallest side of a square of plots that will contain M plots
jac     = mdl.jacobian;
resid   = mdl.residual;
mse     = mean(resid.^2);
rmse    = sqrt(mse);
normJac = jac/rmse;
invCov = normJac'*normJac;
cov = inv(invCov);
% Calculate the covariance matrix from the jacobian. The unconventional
% approach using a single, aveage value for the RMSE and a normalized
% Jacobian is beacuse otherwise the data scales can become unworkable. A
% more conventional and accurate appraoch to calcaulting the covariance
% matrix could be used if all fits were to scaled data. But this isn't
% guaranteed with the code at-present (Oct 2015)
coeff = mdl.x;

f_ = findobj('name', 'Covariance pair plots');
delete(f_);
f_ = figure('name', 'Covariance pair plots');
a_ = zeros(M,1); % handles to plots
ind = 1; % index of current plot
for i = 2:N
    for j = 1:(i-1)
        a_(ind) = subplot(L,L,ind);
        sig_i = sqrt(cov(i,i)); % standard dev of 1st parameter
        sig_j = sqrt(cov(j,j)); % standard dev of 2nd parameter
        nPts = 100;
        x = linspace(coeff(i)-2.5*sig_i, coeff(i)+2.5*sig_i,nPts);
        y = linspace(coeff(j)-2.5*sig_j, coeff(j)+2.5*sig_j,nPts);
        % arrays of values to use for plotting
        P = zeros(nPts,nPts); % matrix of probability density function values
        for u = 1:nPts;
            for v = 1:nPts;
                r = coeff; % vector of parameters at best fit valuves
                r(i) = x(u);
                r(j) = y(v);
                % change the values of the 2 parameters being scanned
                CF = (r-coeff)*invCov*(r-coeff)';
                if CF < 0
                    disp('huh')
                end
                % the covariance form
                P(u,v) = exp(-CF/2)/(2*pi)^(N/2)/sqrt(det(cov));
                % probability density function
            end
        end
        if isreal(P) & sum(sum(isnan(P)))==0
            imagesc(x*PFit(i), y*PFit(j), P);
            set(a_(ind), 'ydir', 'normal');
        else
            text(0.5,0.5, 'No', 'Units', 'Normalized', 'Color', 'r');
        end
        xlabel(PFitNames{i}); ylabel(PFitNames{j});
        ind = ind+1;
    end
end
annotation(f_,'textbox',...
    [0.02 0.94 0.95 0.06],...
    'String',{'Plotting the probability density function using the covariance matrix for the fitted paramters'},...
    'FitBoxToText','off');
