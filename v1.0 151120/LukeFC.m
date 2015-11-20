function n = LukeFC(t,y,pump,P)

% function n = LukeFC(t,y,pump,P)
%
% Function evaluates the excess carrier concentration as a function of time
% and position through the wafer after a light pulse as solved by Luke and
% Cheng, JAP 61, 2282 (1987).
% 
% t = timeseries (nanoseconds)
% y = series of positions (in um) through wafer thickness; front at y=-thick/2,
%   rear at y=thick/2
% pump = selects pump profile: 'Delta', 'Square', or 'Gaussian'
% P = parameters. Parameter set depends on pump profile
%       for Delta:      P = [tau, SRV, thick, alpha, R, difu, N, NaN,   NaN]
%       for Square:     P = [tau, SRV, thick, alpha, R, difu, N, NaN,   T]
%       for Gaussian:   P = [tau, SRV, thick, alpha, R, difu, N, sigma, T]
% tau = bulk lifetime (microseconds)
% SRV = surface recombination velocity (cm/sec)
% thick = wafer thickness (microns)
% alpha = optical absorption (1/cm)
% R = reflectance (scale of 0-1)
% difu = minority carrier diffusion coefficient (cm^2/s)
% N = total pump fluence (#/cm^2)
% T = square pump duration or Gaussian time offset (nanoseconds)
% sigma = Gaussian pump half width (nanoseconds)
%
% n = results array (position, time)
%
% Sept 2014, RJ. Adadpted from LukeAvgFC.m

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

sumSize = 40;
% maximum value of n to terminate the trig series. Can terminate at a small
% value (e.g. 6) for nearly all (y,t) pairs. But for a impulse, the curve
% at t=0 needs many terms in order to converge.
alphaN = zeros(sumSize,1);
betaN = zeros(sumSize,1);
for j = 1:sumSize
    alphaN(j) = fminbnd(@(x) abs(difu*x/SRV-cot(x*thick/2)), 2*pi*(j-1)/thick, 2*pi*j/thick);
    betaN(j) = fminbnd(@(x) abs(difu*x/SRV+tan(x*thick/2)), (2*pi*(j-1)+pi/2)/thick, 2*pi*j/thick);
    % range for fminbnd() to find betaN is a cludge; add pi/2 to avoid
    % finding the local minimum on the left side of the sinularity within
    % each interval
end

% % Plot the roots:
% figure;
% xTmp = linspace(0,2*pi*sumSize,10000);
% subplot(1,2,1)
% plot(xTmp,[cot(xTmp/2); difu*xTmp/thick/SRV]);
% hold all
% plot(alphaN*thick,difu*alphaN/SRV,'ko');
% set(gca,'ylim',[-2*max(abs(difu*alphaN/SRV)), 2*max(abs(difu*alphaN/SRV))]);
% set(gca,'xlim',[min(xTmp) max(xTmp)],'xtick',0:2*pi:max(xTmp))
% grid on
% subplot(1,2,2)
% plot(xTmp,[tan(xTmp/2); -difu*xTmp/thick/SRV]);
% hold all
% plot(betaN*thick,-difu*betaN/SRV,'ko');
% set(gca,'ylim',[-2*max(abs(difu*betaN/SRV)), 2*max(abs(difu*betaN/SRV))]);
% set(gca,'xlim',[min(xTmp) max(xTmp)],'xtick',0:2*pi:max(xTmp))
% grid on

rAlphaN = 1/tau + alphaN.^2*difu;
rBetaN = 1/tau + betaN.^2*difu;
% handy rate expressions
AN = 4*g0p*alphaN*exp(-alpha*thick/2)*(1+R*exp(-alpha*thick))./(alpha^2+alphaN.^2)...
    ./(alphaN*thick+sin(alphaN*thick))...
    .*(alpha*sinh(alpha*thick/2)*cos(alphaN*thick/2)+alphaN*cosh(alpha*thick/2).*sin(alphaN*thick/2));
BN = -4*g0p*betaN*exp(-alpha*thick/2)*(1-R*exp(-alpha*thick))./(alpha^2+betaN.^2)...
    ./(betaN*thick-sin(betaN*thick))...
    .*(alpha*cosh(alpha*thick/2)*sin(betaN*thick/2)-betaN*sinh(alpha*thick/2).*cos(betaN*thick/2));
n = zeros(sumSize,length(y),length(t));
% 3D array to store calculated results. 
%   1st dimension = terms in expansion
%   2nd dimension = thickness
%   3rd dimension = time

y = sort(y);
gInWafer = (y >= -thick/2) & (y <= thick/2);
[~,gIndMin] = min(y(gInWafer));
[~,gIndMax] = max(y(gInWafer));
% indices of start and end positions in the array 'y' that are inside the wafer

switch pump
    case 'Delta'
        tNeg = t < 0;
        for j = 1:sumSize
           for k = gIndMin:gIndMax % step through positions in the wafer
               n(j,k,~tNeg) = AN(j)*cos(alphaN(j)*y(k))*exp(-rAlphaN(j)*t(~tNeg))...
                   +BN(j)*sin(betaN(j)*y(k))*exp(-rBetaN(j)*t(~tNeg));
           end
        end
        g = sum(sum(isnan(n),3),2) == 0;
        g = g & (sum(sum(isinf(n),3),2) == 0);
        % logical array of length sumSize indicating which terms in the sum are
        % neither Inf or NaN
        n = n(g,:,:);
        % only include terms in series that have no NaN or Inf terms
    case 'Square'
        T = P(9); % pulse duration
        g1 = t >= 0 & t <= T; g2 = t > T;
        for j = 1:sumSize
            for k = gIndMin:gIndMax % step through positions in the wafer
                n(j,k,g1) = AN(j)*cos(alphaN(j)*y(k))*(1-exp(-rAlphaN(j)*t(g1)))/rAlphaN(j)/T...
                    +BN(j)*sin(betaN(j)*y(k))*(1-exp(-rBetaN(j)*t(g1)))/rBetaN(j)/T;
                n(j,k,g2) = AN(j)*cos(alphaN(j)*y(k))*(exp(rAlphaN(j)*T)-1)*exp(-rAlphaN(j)*t(g2))/rAlphaN(j)/T...
                    +BN(j)*sin(betaN(j)*y(k))*(exp(rBetaN(j)*T)-1)*exp(-rBetaN(j)*t(g2))/rBetaN(j)/T;
            end
        end
        g = sum(sum(isnan(n),3),2) == 0;
        g = g & (sum(sum(isinf(n),3),2) == 0);
        % logical array of length sumSize indicating which terms in the sum are
        % neither Inf or NaN
        n = n(g,:,:);
        % only include terms in series that have no NaN or Inf terms
    case 'Gaussian'
        sigma = P(8); % pulse half width
        T = P(9); % pulse time offset in nanoseconds
        for j = 1:sumSize
            for k = gIndMin:gIndMax
                n(j,k,:) = AN(j)*cos(alphaN(j)*y(k))*exp((sigma*rAlphaN(j)/sqrt(2))^2)...
                    *exp(-rAlphaN(j)*(t-T)).*erfc((sigma^2*rAlphaN(j)+T-t)/sqrt(2)/sigma)...
                    +BN(j)*sin(betaN(j)*y(k))*exp((sigma*rBetaN(j)/sqrt(2))^2)...
                    *exp(-rBetaN(j)*(t-T)).*erfc((sigma^2*rBetaN(j)+T-t)/sqrt(2)/sigma);
            end
        end
        n = n/2;
        g = sum(sum(isnan(n),3),2) == 0;
        g = g & (sum(sum(isinf(n),3),2) == 0);
        % logical array of length sumSize indicating which terms in the sum are
        % neither Inf or NaN
        n = n(g,:,:);
        % only include terms in series that have no NaN or Inf terms
    otherwise
        n = NaN;
end
n = sum(n,1);
% sum the terms in the expansion
n = reshape(n,length(y),length(t));
% reshape into a 2D array
n = 1e12*n; % convert to cm^-3
end




