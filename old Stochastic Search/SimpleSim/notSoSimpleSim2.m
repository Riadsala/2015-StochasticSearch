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
for beta = 1:3
    beta
    for rcase = 1:3
        results = [];
        r=rvalues(rcase);
        phi = rand*360;
        [x y] = pol2cart(phi, r);
        x = ceil(x); y = ceil(y);
        targetloc = [x y]+N/2;
        for trial=1:numbertrials
            t=0;
            targetfound = 0;
            fix = [];
            fix(1,:) = init_fix;
            trialamps = [];
            while targetfound == 0
                t=t+1;
                tw = min(t, 50);
                dist2target  = PixelsToVisualAngle(norm(fix(t,:)-targetloc));
                probofdetection = linear_model(betaS(beta), dist2target);
                chance = rand*100;
                if chance<probofdetection
                    targetfound = 1;
                    fix(t+1,:) = targetloc;
                    saccamps = [saccamps, dist2target];
                    trialamps = [trialamps, dist2target];
                else
                    % keep drawing random saccades until we get on that falls within the
                    % search area
                    validfix=0;
                    while validfix == 0
                        %                         dist = (wblinv(rand, weibull_a(tw), weibull_b(tw)));
                        dist = find(saccadedata(tw,:)>rand, 1)./2;
                        direction = theta(find(dir_a*rand<dircumv, 1));
                        [x y] = pol2cart(direction, VisualAngleToPixels(dist));
                        saccamps = [saccamps, dist];
                        saccdirs = [saccdirs direction];
                        trialamps = [trialamps, dist];
                        x=round(x)+fix(t,1);
                        y=round(y)+fix(t,2);
                        validfix = (0<x).*(x<=1024).*(0<y).*(y<=1024);
                    end
                    fix(t+1,:) = [x y];                    
                end                
            end            
            results(trial) = t+1;
            trialcoverage = coverage2(fix, 1024);
            for fixn = 1:50;
                if t>(fixn+2)
                    SaccadeTimeSeries{fixn} = [SaccadeTimeSeries{fixn}, trialamps(fixn)];
                end
            end
            
            
        end
        meancoverage(beta, rcase) = mean(trialcoverage);
        meansaccades(beta,rcase) = mean(results(isfinite(results)));
        stdsaccades(beta, rcase) = std(results(isfinite(results)));
    end
end
csvwrite('meansacc.csv', meansaccades);
csvwrite('stdnsacc.csv', stdsaccades);
csvwrite('meancoverage.csv', meancoverage);
figure;
subplot(2,1,1);
hist(saccamps, 0.25:0.5:19.75);
subplot(2,1,2);
rose(saccdirs);
figure
set(gcf, 'DefaultAxesColorOrder', [1 0 0; 0 1 0; 0 0 1]);
errorbar(repmat([1.6,1.65,1.7]',1,3), meansaccades, stdsaccades./sqrt(numbertrials));
legend('r=150', 'r-225', 'r=350', 'location', 'northeast')
axis([1.55 1.75 0 50])


figure;
for fixn=1:50;
    meanSaccAmp(fixn) = mean(SaccadeTimeSeries{fixn});
     stdSaccAmp(fixn) = std(SaccadeTimeSeries{fixn})./size(SaccadeTimeSeries{fixn},2);
     sizeSaccAmp(fixn) = size(SaccadeTimeSeries{fixn},2);
end

errorbar(meanSaccAmp, stdSaccAmp);
xlabel('fixation number (within a trial)')
ylabel('mean saccade amplitude');
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