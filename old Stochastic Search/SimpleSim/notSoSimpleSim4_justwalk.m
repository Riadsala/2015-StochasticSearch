function notSoSimpleSim
close all


% search area dimensions and properties
N = 1024;
init_fix = [N/2, N/2];
numbertrials = 6;
% initialise matrices to store fixation info.
SaccsByPos = cell(5,5);
for trial=1:numbertrials
    fix = [];
    fix(1,:) = init_fix;
    trialamps = [];
    for t=1:30
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
        [d theta] = GetSaccadeDistsAndAngles(x, y) 
        SaccAngles(t) = theta;
        SaccDists(t) = d;        
        SaccsByPos{a,b} = [SaccsByPos{a,b}; [SaccDists(t), SaccAngles(t),t]];
    end 
    figure
    voronoi(fix(:,1), fix(:,2));
   
end % trial end


% 
% figure;
% plot(meanrefix, '-o')
% title('refix per fix')
% 
figure

topleftSaccades = SaccsByPos{5,1};
topleftSaccades(isnan(topleftSaccades(:,1)),:)= [];
subplot(3,3,1); rose(topleftSaccades(:,2))

toprightSaccades = SaccsByPos{1,1};
toprightSaccades(isnan(toprightSaccades(:,1)),:)= [];
subplot(3,3,3); rose(toprightSaccades(:,2))

botleftSaccades = SaccsByPos{5,5};
botleftSaccades(isnan(botleftSaccades(:,1)),:)= [];

botrightSaccades = SaccsByPos{1,5};
botrightSaccades(isnan(botrightSaccades(:,1)),:)= [];

topSaccades = [SaccsByPos{2,1};SaccsByPos{3,1};SaccsByPos{4,1}]

botSaccades =SaccsByPos{2:4,5};
botSaccades(isnan(botSaccades(:,1)),:)= [];

leftSaccades = [ SaccsByPos{5,2}; SaccsByPos{5,3}; SaccsByPos{5,4}];
leftSaccades(isnan(leftSaccades(:,1)),:)= [];

rightSaccades = [ SaccsByPos{1,2}; SaccsByPos{1,3}; SaccsByPos{1,4}];
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

% 
% 
% figure;
% plot(meansaccades, '-o')
% title('num fixes')

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