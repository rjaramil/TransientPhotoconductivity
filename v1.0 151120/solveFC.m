function [nAvg, nProfile, mode_info] = solveFC(t, tDiscrete, y, pump, P, h, solver, injDep, wantProfile);

% function [nAvg, nProfile, mode_info] = solveFC(t, tDiscrete, y, pump, P, h, solver, injDep, wantProfile);
%
% Solves the dynamics of free carrier diffusion and recombination problem
% for modeling transient free-carrier-absorption data. Can call either an
% analytical or a numerical routine.
%
% The analytical routine (LukeAvgFC.m and LukeFC.m) follows Luke & Cheng,
% JAP 61, 2282-2293 (1987). The numerical routine (numFC.m) numerically
% solves the diffusion equation.
%
% t = time points at which to model the average concentration through the
%   wafer (ns)
%
% tDiscrete = timesteps at which to model the spatial profile (ns)
%
% y = space points at which to model the concentration through
%   the wafer (microns)
%
% pump = selects pump profile: 'Delta', 'Square', or 'Gaussian'
%
% P = parameters. Parameter set depends on pump profile. 
%       for Delta:      P = [tau, SRV, thick, alpha, R, difu, N, ~,     ~]
%       for Square:     P = [tau, SRV, thick, alpha, R, difu, N, ~,     T]
%       for Gaussian:   P = [tau, SRV, thick, alpha, R, difu, N, sigma, T]
%
% h = factor controls y-offset due to sample heating. Units are (carrier
%   concentration) / (fluence), or (1/cm).
%
% solver = selects the solver: 'Analytical', or 'Numerical'
%
% injDep = structure to describe the injection-dependence of the model.
%   Only applies to Numerical solver.
%   injDep.tauModel = selects bulk lifetime model
%   injDep.difuModel = selects diffusion model
%   injDep.effMassRatio = majority/minority carrier effective mass ratio
%   injDep.majConc = majority carrier concentration (cm^-3)
%
% wantProfile = If true, then the carrier concentration profile in space on
%   the grid 'y' is calculated for time points tDiscrete. If false, this data
%   is not returned.
%
% nAvg = average excess free carrier concentration through the wafer as a
%   function of time points 't' (cm^-3)
%
% nProfile = spatial profile of excess free carrier concentration through
%   the water at time points 'tDiscrete' and spatial points 'y' (cm^-3).
%   Size(nProfile) = (m, n) where m = length(y) and n = length(tDiscrete)
%
% mode_info = details of the harmonic modes for the analytical solution.
%
% July 2015, R Jaramillo

switch solver
    
    case 'Analytical'
        
        [nAvg, mode_info] = LukeAvgFC(t, pump, P);
        if wantProfile
            nProfile = LukeFC(tDiscrete, y, pump, P);
        else
            nProfile = NaN;
        end
        
    case 'Numerical'
        
        [n, nAvg] = numFC(t, y, pump, P, injDep);
        mode_info = NaN;
        if wantProfile
            nProfile = interp2(t, y', n, tDiscrete, y');
            % interpolate full solution to find profile at only times
            % 'tDiscrete'
        else
            nProfile = NaN;
        end
        
    otherwise
        nAvg = NaN; nProfile = NaN; mode_info = NaN;
        
end

% Add the heating contribution:

N = P(7); 
% N = total pump fluence (#/cm^2)
nOffset = zeros(size(nAvg));
switch pump
    
    case 'Delta'
        
        tNeg = t < 0;
        nOffset(~tNeg) = h*N;       
        
    case 'Square'
        
        T = P(9);
        % T = square pump duration or Gaussian time offset (nanoseconds)
        g1 = t >= 0 & t <= T; g2 = t > T;
        nOffset(g1) = h*t(g1)*N/T;
        nOffset(g2) = h*N;
    
    case 'Gaussian'
        
        T = P(9);
        % T = square pump duration or Gaussian time offset (nanoseconds)
        sigma = P(8);
        % sigma = Gaussian pump half width (nanoseconds)
        nOffset(:) = (h*N/2)*(1 + erf((t-T)/sqrt(2)/sigma));

end
nAvg = nAvg + nOffset;

if (isrow(t) & ~isrow(nAvg)) | (iscolumn(t) & ~iscolumn(nAvg))
    nAvg = nAvg';
end
% makes sure that input and output are the same size

end
