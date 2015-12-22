function notSoSimpleSim
close all

% search area dimensions and properties
N = 1024;
betaS = [1.6, 1.65, 1.7];
rvalues = [100 225 350];
init_fix = [N/2, N/2];
numbertrials = 100;
% initialise matrices to store fixation info.
SaccsByPos = cell(5,5);
saccamps = [];
saccdirs = [];
SaccadeTimeSeries = cell(1,50);
for beta =1:3
    %comment this to work over all r
    results = [];
    for rcase = 1:3
        %uncomment this to work over all r
        %         results = zeros(1,numbertrials);
        for trial=1:numbertrials
            r=rvalues(rcase);
            phi = rand*360;
            [x y] = pol2cart(phi, r);
            x = ceil(x); y = ceil(y);
            targetloc = [x y]+N/2;
            t=0;
            targetfound = 0;
            fix = [];
            fix(1,:) = init_fix;
            trialamps = [];
            while targetfound == 0
                t=t+1;
                dist2target  = PixelsToVisualAngle(norm(fix(t,:)-targetloc));
                probofdetection = linear_model(betaS(beta), dist2target);
                chance = rand*100;
                if chance<probofdetection
                    targetfound = 1;
                    fix(t+1,:) = targetloc; %#ok<AGROW>
                else
                    % keep drawing random saccades until we get on that falls within the
                    % search area
                    validfix=0;
                    while validfix == 0
                        % work out if we're in a corner, edge of middle
                        % location
                        [l a b] = getBin(fix(t,1),fix(t,2));
                        [x y] = getSaccadeFromDist3(l,t);
                        xf=x+fix(t,1);
                        yf=y+fix(t,2);
                        validfix = (0<xf).*(xf<=1024).*(0<yf).*(yf<=1024);
                    end
                    fix(t+1,:) = [xf yf]; %#ok<AGROW>
                    [d theta] = GetSaccadeDistsAndAngles(x, y);
                    SaccAngles(t) = theta;
                    SaccDists(t) = d;
                    SaccsByPos{a,b} = [SaccsByPos{a,b}; [SaccDists(t), SaccAngles(t),t]];
                end
            end % trial end            
%             close all
%             figure
%             plot(fix(:,1), fix(:,2), 'b-');
%             hold on
%             plot(targetloc(1), targetloc(2), 'rx');
%             axis([1 1024, 1 1024]);
%             filename = ['scanpath_b=' num2str(beta) '_r=' int2str(rcase) ...
%                 '_trial=' int2str(trial) '.jpg'];
%             print('-djpeg', filename);
            results(trial) = t+1;
            numrefix(trial) = checkForReFixations(fix);
        end % trial end
        meansaccades(beta,rcase) = mean(results(isfinite(results)));
        stdsaccades(beta,rcase) = std(results(isfinite(results)));
        meanrefix(beta,rcase) = mean(numrefix./results);
        
    end
    
end
csvwrite('meansacc.csv', meansaccades);
csvwrite('stdnsacc.csv', stdsaccades);
mean(meanrefix,2)

 csvwrite('meanrefix.csv', meanrefix);
figure
plot(meansaccades)
figure;
plot(meanrefix, '-o')
title('refix per fix')

% figure('name', 'number of fixations by time and position')
% plotctr=0;
% for b=1:5
%     for a=5:-1:1
%         clear S fixnums
%         plotctr = plotctr+1;
%         subplot(5,5, plotctr)
%         S= SaccsByPos{a,b};
%         S(isnan(S(:,1)),:)= [];
%         for t=1:30
%             fixnums(t) = size(find(S(:,3)==t),1);
%             SaccAmpsByTimeAndPos(a,b,t) = mean(S(S(:,3)==t),1);
%         end
%         plot(fixnums);
%         axis([0 30 0 40])
%     end
% end
% figure('name', 'effect of time on saccade amplitude');
% plotctr=0;
% for b=1:5
%     for a=5:-1:1
%         clear T tmp tmpnan
%         plotctr = plotctr+1;
%         subplot(5,5, plotctr)
%         tmp =reshape(SaccAmpsByTimeAndPos(a,b,:), [1,30]);
%         plot(tmp);
%         tmpnan=isnan(tmp);
%         T = 1:30;
%         T(tmpnan)=[];
%         tmp(tmpnan) = [];
%         regmodels(a,b,:) = polyfit(T,tmp,1);
%         axis([0 30 0 15])
%     end
% end


%%%%%%%%%
% plot mean sacc by position (not time)
%%%%%%%%%
% figure('name', 'mean saccade amplitude by position')
%
% for b=1:5
%     for a=5:-1:1
%         clear S fixnums
%         S= SaccsByPos{a,b};
%         S(isnan(S(:,1)),:)= [];
%         meansacc(a,b) = mean(S(:,1));
%     end
% end
% bar3(meansacc)

topleftSaccades = SaccsByPos{5,1};
topleftSaccades(isnan(topleftSaccades(:,1)),:)= [];

toprightSaccades = SaccsByPos{1,1};
toprightSaccades(isnan(toprightSaccades(:,1)),:)= [];

botleftSaccades = SaccsByPos{5,5};
botleftSaccades(isnan(botleftSaccades(:,1)),:)= [];

botrightSaccades = SaccsByPos{1,5};
botrightSaccades(isnan(botrightSaccades(:,1)),:)= [];

topSaccades = SaccsByPos{2:4,1};
topSaccades(isnan(topSaccades(:,1)),:)= [];

botSaccades =SaccsByPos{2:4,5};
botSaccades(isnan(botSaccades(:,1)),:)= [];

leftSaccades = SaccsByPos{5,2:4};
leftSaccades(isnan(leftSaccades(:,1)),:)= [];

rightSaccades = SaccsByPos{1,2:4};
rightSaccades(isnan(rightSaccades(:,1)),:)= [];

midSaccades = [SaccsByPos{2,2}; SaccsByPos{3,2}; SaccsByPos{4,2};...
    SaccsByPos{2,3}; SaccsByPos{3,3}; SaccsByPos{4,3};...
    SaccsByPos{2,4}; SaccsByPos{3,4}; SaccsByPos{4,4}];
midSaccades(isnan(midSaccades(:,1)),:)= [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flip corners around to get composite
figure
CornerTL = [topleftSaccades; fliphorizontally(toprightSaccades);
    flipvertically(botleftSaccades); flipvertically(fliphorizontally(botrightSaccades))];
subplot(2,2,1);
rose(CornerTL(:,2));
title('Corner Behaviour');
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mirror horizontal edge behaviour in vertical axis, then flip bottom edge and merge
% with top
EdgeT = [topSaccades;fliphorizontally(topSaccades);...
    flipvertically(botSaccades); flipvertically(fliphorizontally(botSaccades))];
subplot(2,2,2)
rose(EdgeT(:,2));
title('Horizontal Edge')
% and similar for the vertical edge
EdgeL = [leftSaccades; flipvertically(leftSaccades);...
    fliphorizontally(rightSaccades); fliphorizontally(flipvertically(rightSaccades))];
subplot(2,2,3)
rose(EdgeL(:,2));
title('Vertical Edge')
%and display middle
subplot(2,2,4);
rose(midSaccades(:,2));
title('Central Region')

    function S=fliphorizontally(S)
        S(:,2) = -(S(:,2)-pi/2)+pi/2;
        S(:,2) = mod(S(:,2),2*pi);
    end

    function S=flipvertically(S)
        S(:,2) = -(S(:,2)-pi/4)+pi/4;
        S(:,2) = S(:,2)+3*pi/2;
        S(:,2) = mod(S(:,2),2*pi);
    end

figure
subplot(2,2,1); map = SurfaceFromPoints(CornerTL);

subplot(2,2,2); map = SurfaceFromPoints(EdgeT);

subplot(2,2,3); map = SurfaceFromPoints(EdgeL);

subplot(2,2,4); map = SurfaceFromPoints(midSaccades);



figure;
plot(meansaccades, '-o')
title('num fixes')

figure('Position',[1 1 1000 1000])
plotctr=0;
for b=1:5
    for a=5:-1:1
        clear S
        plotctr = plotctr+1;
        subplot(5,5, plotctr)
        S= SaccsByPos{a,b};
        [Sx Sy] = pol2cart(S(:,2), S(:,1));
        S = [Sx Sy];
        S(isnan(S(:,1)),:)= [];
        size(S)
        ContourPlotFromPoints(S)
    end
end


   function [d theta] = GetSaccadeDistsAndAngles(xdiff, ydiff)        
        d = PixelsToVisualAngle(sqrt(xdiff^2+ydiff^2));
        if xdiff == 0
            theta = pi/2;
        elseif xdiff>0
            theta = atan(ydiff/xdiff);
        elseif xdiff <0
            theta = pi+atan(ydiff/xdiff);
        end        
    end % SaccadeInfo
end


function [l, a,b] = getBin(x,y)
a = 0;
b = 0;
if (x>0)&&(x<205)
    a=1;
elseif (x>=205)&&(x<410)
    a=2;
elseif (x>=410)&&(x<615)
    a=3;
elseif (x>=615)&&(x<820)
    a=4;
elseif (x>=820)&&(x<1025)
    a=5;
end
if (y>0)&&(y<205)
    b=1;
elseif (y>=205)&&(y<410)
    b=2;
elseif (y>=410)&&(y<615)
    b=3;
elseif (y>=615)&&(y<820)
    b=4;
elseif (y>=820)&&(y<1025)
    b=5;
end
if (a==1)&&(b==1)
    l=[1 1];
elseif (a==1)&&(b==5)
    l=[1 2];
elseif (a==5)&&(b==1)
    l=[1 3];
elseif (a==5)&&(b==5)
    l=[1 4];
elseif ((b==2)+(b==3)+(b==3))&&(a==1)
    l=[2 1];
elseif ((b==2)+(b==3)+(b==3))&&(a==5)
    l=[2 2];
elseif ((a==2)+(a==3)+(a==3))&&(b==1)
    l=[3 1];
elseif ((a==2)+(a==3)+(a==3))&&(b==5)
    l=[3 2];
else
    l = [4 1];
end
end



function p = linear_model(beta, r)
p = -597 + 409.145*beta - 11.441*r;
% switch beta
%     case 1.6
%         p = 42.8.*exp(r.*(-.274));
%     case 1.65
%       p = 132.9.*exp(r.*(-.41));
%     case 1.7
%         p = 150.9.*exp(r.*(-.317));
% end


end