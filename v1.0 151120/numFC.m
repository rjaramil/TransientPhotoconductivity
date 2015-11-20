function [nOut, nAvgOut] = numFC(t, y, pump, P, injDep)

% function [nOut, nAvgOut] = numFC(t, y, pump, P, injDep)
%
% Function evaluates the thickness-averaged excess minority carrier
% concentration as a function of time after a light pulse. The average is
% through the thickness of the wafer. The solution is calculated
% numerically. Unlike the analytic solution of Luke & Cheng (1987), this
% allows for nonlinear terms (e.g. injection dependence)
%
% t = time points (ns) at which to solve the model
%
% y = spatial points (microns) at which to solve the model
%
% pump = selects pump profile: 'Delta', 'Square', or 'Gaussian'
%
% P = parameters. Parameter set depends on pump profile
%       for Delta:      P = [tau, SRV, thick, alpha, R, difu, N, NaN,   NaN]
%       for Square:     P = [tau, SRV, thick, alpha, R, difu, N, NaN,   T]
%       for Gaussian:   P = [tau, SRV, thick, alpha, R, difu, N, sigma, T]
%
% injDep = structure to describe the injection-dependence of the model.
%   Only applies to Numerical solver.
%   injDep.tauModel = selects bulk lifetime model
%   injDep.difuModel = selects diffusion model
%   injDep.effMassRatio = majority/minority carrier effective mass ratio
%   injDep.majConc = majority carrier concentration (cm^-3)
%
% nOut = solution as a function of time and space. Size(nOut) = m x n where
% m = length(y) and n = length(t)
%
% nAvgOut = thickness-averaged solution as a function of time
%
% July 2015, R. Jaramillo


%%
% Read in parameters. Make all lengths in microns, all times in
% nanoseconds
tau = P(1); % bulk lifetime (nanoseconds)
SRV = 1e-5*P(2); % SRV (microns/nanosecond)
thick = P(3); % thickness (microns)
alpha = 1e-4*P(4); % absorption coefficient (1/microns)
R = P(5);
difu = 1e-1*P(6); % diffusion (microns^2/nanoseconds)
N = 1e-8*P(7); % fluence (1/microns^2)

g0pp = N*alpha*(1-R);
g0p = g0pp/(1-(R*exp(-alpha*thick))^2);
g0 = N*alpha*(1-R)/(1-R*exp(-alpha*thick));
% g_0: an average generation

%%
% Establish the computational grid in normalized units:

numYPts = min(floor(thick/0.005), floor(thick*30/(1/alpha)));
% % approximate stepsize: 30 steps per absorption length or 5 nm steps,
% % whichever is larger
dy = thick/(numYPts - 3);
yPts = linspace((-thick/2 -dy), (thick/2 + dy), numYPts);
% gridpoints in space

% Establish start and stop times and timestep. Make sure that start and
% stop times span the generation pulse, and that the timestep is small
% enough to represent the pulse:

switch pump
    case 'Delta'
        tStart = 0;
    otherwise
        tStart  = min(t);
end
% for delta function pump, pump occurs at time=0 and the starting time
% carrier profile is explicitly entered
tStop   = max(t);
dt = 100e-6; % start with 100 fs timesteps
numTPts = floor((tStop - tStart)/dt) + 1;
% number of time points
b = difu*dt/(dy)^2;
% stability parameter. Using the largest possible value of the diffusivity,
% assumed here to be the minority carrier (electron) diffusivity. This
% needs to be generalized for different semiconductors.
while b > 0.3
    numTPts = numTPts*2;
    dt = (tStop - tStart)/(numTPts-1);
    b = difu*dt/(dy)^2;
end
% double the number of timesteps until the stability criterion is satisfied
tPts = linspace(tStart, tStop, numTPts);
% gridpoints in time


%%
% Create solution array, and initialize the values at the first position in
% time

n = zeros(numYPts, numTPts);
% array to store all of the solutions

% In case of delta pump, initialize values. Otherwise they all start at
% zero.
if strcmp(pump, 'Delta')
    n(:, 1) = g0pp*(exp(-alpha*(yPts + thick/2)) ...
                + P(5)*exp(-alpha*thick)*exp(-alpha*(-yPts + thick/2)))...
                /(1 - (P(5)*exp(-alpha*thick))^2);
end
    

%%
% Step forward in time
wb_ = waitbar(0, 'Calculating...', 'name', 'numFC.m progress');
for k = 2:numTPts
    
    waitbar(k/numTPts)

    % Build the problem

    M = zeros(numYPts, numYPts);
    % the problem matrix
    
    % Populate the problem matrix for all but the 1st and last rows:
    for u = 2:(numYPts - 1)
        for v = 1:numYPts
            if v == (u - 1)
                M(u, v) = difu_n(n(u, k-1))*dt/(dy)^2;
            elseif v == u
                M(u, v) = 1 - 2*difu_n(n(u, k-1))*dt/(dy)^2 - dt/tau_n(n(u, k-1));
            elseif v == (u + 1)
                M(u, v) = difu_n(n(u, k-1))*dt/(dy)^2;
            end
        end
    end
    
    g = zeros(numYPts, 1);
    % the generation matrix
    
    % Populate the generation matrix:
    for u = 1:numYPts
       g(u) = dt*Gen(u, k-1);
    end
    
    n(:, k) = M*n(:, k-1) + g;
    % calculate all but the ghost points
    
    n(1, k)     = n(3, k) - (2*dy*SRV_n(n(2, k))/difu_n(n(2, k)))*n(2, k);
    n(end, k)   = n(end-2, k) - (2*dy*SRV_n(n(end-1, k))/difu_n(n(end-1, k)))*n(end-1, k);
    % calculate the ghost points based on time=k values

end
delete(wb_)

n = 1e12*n;
% convert to 1/cm^3

nAvg = mean(n(2:end-1, :), 1);
% average through the sample. Calculate this average first, before
% interpolating the spatial grid, just in case the input spatial grid
% didn't cover the full sample thickness

switch pump

    case 'Delta'
        nAvgOut = zeros(size(t)); % array of output values
        g = t >= 0;
        nAvgOut(g) = interp1(tPts, nAvg, t(g));

    otherwise
        nAvgOut = interp1(tPts, nAvg, t);

end
% put the average values on the input time grid. In case of Delta function
% pump, explicitly set the values to 0 for t<0

nOut = interp2(tPts, yPts', n, t, y');
% interpolate to find values of solution at the input time and space points

%%
% Output figures

% % Figure comparing numerical solution with analytic solution of Luke and
% % Cheng
% figure;
% subplot(1,2,1)
% plot(t, nAvgOut, '.-', 'displayname', 'numFC')
% xlabel('time (ns)')
% ylabel('average n (cm^{-3})')
% [nLukeAvg, ~] = LukeAvgFC(t, pump, P);
% hold all
% plot(t, nLukeAvg, '.-', 'displayname', 'LukeFC')
% title(sprintf('dt = %.2e (ns), dy = %.2e \\mum', dt, dy));
% subplot(1,2,2)
% plot(t, nAvgOut./nLukeAvg, '.-')
% title('numerical avg / Luke avg')
% xlabel('time (ns)')


%%
% The function tau_n(nIn) defines how the bulk lifetime varies with injection
% nIn = minority carrier concentration

function out = tau_n(nIn)
    
    switch injDep.tauModel
        
        case 'Constant'
            
            out = tau;
            
        case 'Simple SRH'
            out = tau*(1 + 2*1e12*nIn/injDep.majConc)/(1 + 1e12*nIn/injDep.majConc);
            % lifetime depends on injection level. Careful to convert value
            % 'nIn' to units of cm^-3
            
        otherwise
            
            out = tau;
            
    end
    
end

%%
% The function difu_n(nIn) defines how the diffusion constant varies with
% injection.
% nIn = minority carrier concentration
% q = 1;
% kB = 8.617e-5; % eV/K
% handles.mu = handles.difu*q/kB/handles.temp;

function out = difu_n(nIn)
    
    switch injDep.difuModel
        
        case 'Constant'
            
            out = difu;
            
        case 'Ambipolar'
            
            out = difu*(injDep.majConc + 2*nIn)/...
                (injDep.majConc + nIn*(1 + injDep.effMassRatio));
            
        otherwise
            
            out = difu;
            
    end
end

%%
% The function SRV_n(nIn) defines how the SRV varies with
% injection.
% nIn = minority carrier concentration

function out = SRV_n(nIn)
    out = SRV;
end

%% 
% The function Gen(u, k) defines the carrier generation.
% k = time index
% u = space index

function out = Gen(u, k)
    
    switch pump
        
        case 'Delta'
            out = 0;
        
        case 'Square'
            if tPts(k) < 0
                out = 0;
                
            elseif tPts(k) < P(9)
                out = 1/P(9);
                % time dependence
                out = out*g0pp*(exp(-alpha*(yPts(u) + thick/2)) ...
                + P(5)*exp(-alpha*thick)*exp(-alpha*(-yPts(u) + thick/2)))...
                /(1 - (P(5)*exp(-alpha*thick))^2);
                % space dependence
                
            else
                out = 0;
                
            end
            
        case 'Gaussian'
            out = (1/sqrt(2*pi)/P(8))*exp(-(tPts(k) - P(9))^2/2/P(8)^2);
            % time dependence
            out = out*g0pp*(exp(-alpha*(yPts(u) + thick/2)) ...
                + P(5)*exp(-alpha*thick)*exp(-alpha*(-yPts(u) + thick/2)))...
                /(1 - (P(5)*exp(-alpha*thick))^2);
            % space dependence
            
        otherwise
            out = 0;
    
    end

end

end




