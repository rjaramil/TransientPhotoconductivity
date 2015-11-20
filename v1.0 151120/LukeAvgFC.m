function [n, mode_info] = LukeAvgFC(t,pump,P)

% function [n, mode_info] = LukeAvgFC(t,pump,P)
%
% Function evaluates the thickness-averaged excess minority carrier
% concentration as a function of time after a light pulse as solved by Luke
% and Cheng, JAP 61, 2282 (1987). The average is through the thickness of
% the wafer.
%
% Optionally, could compute the position-dependent excess carrier
% concentration as a function of time, if a series of points is passed in
% varargin.
% 
% t = timeseries (nanoseconds)
% pump = selects pump profile: 'Delta', 'Square', or 'Gaussian'
% P = parameters. Parameter set depends on pump profile
%       for Delta:      P = [tau, SRV, thick, alpha, R, difu, N, ~,     ~]
%       for Square:     P = [tau, SRV, thick, alpha, R, difu, N, ~,     T]
%       for Gaussian:   P = [tau, SRV, thick, alpha, R, difu, N, sigma, T]
% tau = bulk lifetime (nanoseconds)
% SRV = surface recombination velocity (cm/sec)
% thick = wafer thickness (microns)
% alpha = optical absorption (1/cm)
% R = reflectance (scale of 0-1)
% difu = minority carrier diffusion coefficient (cm^2/s)
% N = total pump fluence (#/cm^2)
% T = square pump duration or Gaussian time offset (nanoseconds)
% sigma = Gaussian pump half width (nanoseconds)
%
% n = thickness-averaged excess minority carrier concentration (cm^-3)
%
% April 2014, RJ
%
% August 2014, RJ: Change Gaussian calculation to allow for pump pulse
% centered at any time, including t < 0. See calculations on 8/19/14.
%
% Dec 2014, RJ: Started adding a pump-dependent y-offset to the model to
% account for sample heating.
%
% June 2015, RJ: Added calculation of mode info, returned to the calling
% function
%
% July 2015, RJ: Removed pump-dependent y-offset. This can be added
% outside this function

s = size(t);
if iscolumn(t)
    t = t';
end
% Makes sure that 't' is a row vector

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

sumSize = 10;
% maximum value of n to terminate the trig series
alphaN = zeros(sumSize,1);
for j = 1:sumSize
    alphaN(j) = fminbnd(@(x) abs(difu*x/SRV-cot(x*thick/2)), 2*pi*(j-1)/thick, 2*pi*j/thick);
end
rAlphaN = 1/tau + alphaN.^2*difu;
% handy rate expression
tauAlphaN = 1./rAlphaN;
% equivalent time constant

% % Plot to show the roots:
% x = linspace(0,max(alphaN),1000);
% figure;
% plot(x*thick/2/pi,x*difu/SRV,'displayname','x*D/S')
% hold all
% plot(x*thick/2/pi,cot(x*thick/2),'displayname','cot(x*d/2)')
% plot(alphaN*thick/2/pi,alphaN*difu/SRV,'o','displayname','roots')
% set(gca,'ylim',(difu/SRV)*max(alphaN)*[-1 1]);
% title('roots \alpha_N')
% xlabel('x*d/2/pi')
% legend('show')
% grid on

n = zeros(sumSize,length(t));
mode_amps = zeros(sumSize, 1);

switch pump
    
    case 'Delta'
        tNeg = t < 0;
        for j = 1:sumSize
            amp = (sin(alphaN(j)*thick/2)/(alpha^2+alphaN(j)^2)/(alphaN(j)*thick+sin(alphaN(j)*thick)))...
                *(alpha*sinh(alpha*thick/2)*cos(alphaN(j)*thick/2)+alphaN(j)*cosh(alpha*thick/2)*sin(alphaN(j)*thick/2));
            amp = amp*8*g0*exp(-alpha*thick/2)/thick;
            mode_amps(j) = amp;
            % amplitude for the jth mode
            n(j,~tNeg) = amp*exp(-rAlphaN(j)*t(~tNeg));
        end
        g = sum(isnan(n),2) == 0;
        g = g & sum(isinf(n),2) == 0;
        % logical of size [sumSize,0] indicating which terms in the sum are
        % neither Inf or NaN
        n = n(g, :);
        % only include terms in series that have no NaN or Inf terms

    case 'Square'
        T = P(9); % pulse duration
        g1 = t >= 0 & t <= T; g2 = t > T;
        for j = 1:sumSize
            ampDuring = (sin(alphaN(j)*thick/2)/rAlphaN(j)/(alpha^2+alphaN(j)^2)...
                /(alphaN(j)*thick+sin(alphaN(j)*thick)))*(alpha*sinh(alpha*thick/2)*cos(alphaN(j)*thick/2)...
                +alphaN(j)*cosh(alpha*thick/2)*sin(alphaN(j)*2/thick));
            ampDuring = ampDuring*8*g0*exp(-alpha*thick/2)/thick/T;
            n(j,g1) = ampDuring*(1-exp(-rAlphaN(j)*t(g1)));
            ampAfter = (sin(alphaN(j)*thick/2)*(exp(rAlphaN(j)*T)-1)/rAlphaN(j)/(alpha^2+alphaN(j)^2)...
                /(alphaN(j)*thick+sin(alphaN(j)*thick)))*(alpha*sinh(alpha*thick/2)*cos(alphaN(j)*thick/2)...
                +alphaN(j)*cosh(alpha*thick/2)*sin(alphaN(j)*2/thick));
            ampAfter = ampAfter*8*g0*exp(-alpha*thick/2)/thick/T;
            mode_amps(j) = ampAfter;
            n(j,g2) = ampAfter*exp(-rAlphaN(j)*t(g2));
        end
        g = sum(isnan(n),2) == 0;
        g = g & sum(isinf(n),2) == 0;
        % logical of size [sumSize,0] indicating which terms in the sum are
        % neither Inf or NaN
        n = n(g,:);
        % only include terms in series that have no NaN or Inf terms
        
    case 'Gaussian'
        sigma = P(8); % pulse half width
        T = P(9); % pulse time offset in nanoseconds
        for j = 1:sumSize
            amp = (sin(alphaN(j)*thick/2)/(alpha^2+alphaN(j)^2)/(alphaN(j)*thick+sin(alphaN(j)*thick)))...
                *(alpha*sinh(alpha*thick/2)*cos(alphaN(j)*thick/2)...
                +alphaN(j)*cosh(alpha*thick/2)*sin(alphaN(j)*thick/2))...
                *exp((sigma*rAlphaN(j)/sqrt(2))^2);
            amp = amp*8*g0*exp(-alpha*thick/2)/thick;
            amp = amp/sqrt(2*pi)/sigma;
            amp = amp*sqrt(pi)/2;
            amp = amp*sqrt(2)*sigma;
            mode_amps(j) = amp;
            % amplitude of this spatial mode
            n(j,:) = amp*exp(-rAlphaN(j)*(t-T)).*erfc((sigma^2*rAlphaN(j)+T-t)/sqrt(2)/sigma);
            % multiply by the time dependence
        end
        g = (sum(isnan(n),2) == 0) & (sum(isinf(n),2) == 0);
        % logical of size [sumSize,0] indicating which terms in the sum are
        % neither Inf or NaN
        n = n(g,:);
        % only include terms in series that have no NaN or Inf terms
        
end

n           = 1e12*n;
mode_amps   = 1e12*mode_amps;
% convert to cm^-3


% figure;
% axes
% hold all
% ar_ = area(t, n');
% xlabel('time (ns)')
% ylabel('cm^{-3}')
% title('area plot of modes')

n = sum(n,1);
% sum the terms in the expansion

mode_info = {tauAlphaN mode_amps};
% mode info to return

% figure;
% subplot(2,1,1)
% plot(1:sumSize, tauAlphaN, 'o-')
% set(gca,'xtick',1:sumSize,'xlim',[1 sumSize])
% xlabel('mode #')
% ylabel('mode time constant (ns)')
% subplot(2,1,2)
% plot(tauAlphaN, mode_amps, 'o-')
% xlabel('mode time constant (ns)')
% ylabel('mode amplitude (1/cm^3')

end




