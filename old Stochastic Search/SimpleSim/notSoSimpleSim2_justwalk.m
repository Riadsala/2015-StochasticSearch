function notSoSimpleSim
close all

rand('state', 1);


% get Directional Cumulative Probability Distribution from Empirical Data
thetares = 0.01;
[dircumv, dir_a] = getDirCumv(thetares);
theta = 0:thetares:2*pi;

% Saccade Amplitude uses raw Empirical Data
saccadedata = csvread('SaccadeAmpFixInfo.csv');

% search area dimensions
N = 1024;


init_fix = [N/2, N/2];

% initialise matrices to store fixation info.


betaS = [1.6, 1.65, 1.7];
rvalues = [100 225 350];
numbertrials = 100;


saccamps = [];
saccdirs = [];
SaccadeTimeSeries = cell(1,50);
V_Area = [];
for trial = 1:105
    fix = zeros(30,2);
for t=1:30      
    % keep drawing random saccades until we get on that falls within the
    % search area
    validfix=0;
    while validfix == 0
        %                         dist = (wblinv(rand, weibull_a(tw), weibull_b(tw)));
        dist = find(saccadedata(t,:)>rand, 1)./2;
        direction = theta(find(dir_a*rand<dircumv, 1));
        [x y] = pol2cart(direction, VisualAngleToPixels(dist));
        saccamps = [saccamps, dist];
        saccdirs = [saccdirs direction];
%         trialamps = [trialamps, dist];
        x=round(x)+fix(t,1);
        y=round(y)+fix(t,2);
        validfix = (0<x).*(x<=1024).*(0<y).*(y<=1024);
    end
    fix(t+1,:) = [x y];
end
V_Area = [V_Area,v_area(fix)];
end

h=hist(log(V_Area), 100);
h = h./size(V_Area,2);

bar(1:1.2:120,h, 1)
xlabel('log(Area of Voronoi Cell)');
axis([0 120 0 0.09]);






end


function [dircumv, dir_a] = getDirCumv(thetares)
P = csvread('FourierComponents.csv');
theta = 0:thetares:2*pi;
dircumv = P(1)*theta;
for j=2:64;
    dircumv =  dircumv + (-1)^(j-1).*real(P(j)).*sin((j-1)*theta)./(j-1) + (-1)^j*imag(P(j)).*cos((j-1)*theta)./(j-1);
end
dircumv = dircumv./64;
dir_a = max(dircumv);
end

function p = linear_model(beta, r)


p = -597 + 409.145*beta - 11.441*r;

end