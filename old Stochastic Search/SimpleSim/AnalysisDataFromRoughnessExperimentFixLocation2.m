function AnalysisDataFromRoughnessExperimentFixLocation
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% .m file to read, process, and display results from Eye tracking
% experiment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment Information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numberofParts = 4;
rSet = [100, 225, 350];
RMSSet = [0.9, 1.1];
betaSet = [1.6, 1.65, 1.7];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Participants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
People = {'mk', 'my', 'lm', 'IHH', 'hw', 'am', 'ED'};
numberOfPeople = size(People, 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Empty Arrays to Save Results Into
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PData = creatEmptyResultsStructure;
SaccadeDists = [];
SaccadeAngles = [];
IndvSaccadeAmps = cell(50,1);
IndvSaccadeAmpsByBeta = cell(50,3,2);
IndvSaccadeAmps20 = cell(20,1);


% make a 4x4 cell array for sorting saccades by fixation location
SaccsByPos = cell(5,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Each Person, Extract Results from Clearview Export Text Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Pctr = 1:numberOfPeople
    disp(['Now extracting info for person number ' int2str(Pctr)])
    Person = People{Pctr};
    [PData(Pctr).Results] = SortOutEventData(Person, numberofParts);
    for run = 1:numberofParts
        Fixs = getFixesFromTextFile;
        for trialctr = 1:90
            trialnmbr = trialctr+(run-1)*90;
            % for each trial, get fix info
            StartTime = PData(Pctr).Results(trialnmbr, 1);
            EndTime = PData(Pctr).Results(trialnmbr,2);
            Fixations = getFixInfo(Fixs);
            PData(Pctr).Results(trialnmbr,9) = Fixations.number;
            [SaccDists, SaccAngles] = GetSaccadeDistsAndAngles;
            SaccadeAngles = [SaccadeAngles, SaccAngles]; %#ok<AGROW>
            SaccadeDists = [SaccadeDists, SaccDists]; %#ok<AGROW>
            % go through each saccade and bin it by fixation location
            for f=1:(Fixations.number-2) %#ok<FXUP>
                [a b] = getBin(Fixations.x(f), Fixations.y(f));
                if (a>0)&&(b>0)
                    SaccsByPos{a,b} = [SaccsByPos{a,b}; [SaccDists(f), SaccAngles(f),f]];
                end
            end
        end % end trialctr loop
    end % end run loop
end % end Pctr loop

%%%%%%%%%%%
% now plot saccades stats
%%%%%%%%%%%%%%
figure('Position',[1 1 1000 1000])
plotctr=0;
for b=1:5
    for a=5:-1:1
        clear S
        plotctr = plotctr+1;
        subplot(5,5, plotctr)
        S= SaccsByPos{a,b};
        %         [Sx Sy] = pol2cart(S(:,2), S(:,1));
        %         S = [Sx Sy];
        S(isnan(S(:,1)),:)= [];
        ContourPlotFromPoints2(S)
    end
end

    function plotSaccsByTimeAndPosition
        figure('name', 'number of fixations by time and position')
        plotctr=0;
        for b=1:5
            for a=5:-1:1
                clear S fixnums
                plotctr = plotctr+1;
                subplot(5,5, plotctr)
                S= SaccsByPos{a,b};
                S(isnan(S(:,1)),:)= [];
                for t=1:30
                    fixnums(t) = size(find(S(:,3)==t),1);
                    SaccAmpsByTimeAndPos(a,b,t) = mean(S(S(:,3)==t),1);
                end
                plot(fixnums);
                axis([0 30 0 40])
            end
        end
        figure('name', 'effect of time on saccade amplitude');
        plotctr=0;
        for b=1:5
            for a=5:-1:1
                clear T tmp tmpnan
                plotctr = plotctr+1;
                subplot(5,5, plotctr)
                tmp =reshape(SaccAmpsByTimeAndPos(a,b,:), [1,30]);
                plot(tmp);
                tmpnan=isnan(tmp);
                T = 1:30;
                T(tmpnan)=[];
                tmp(tmpnan) = [];
                regmodels(a,b,:) = polyfit(T,tmp,1);
                axis([0 30 0 15])
            end
        end
    end

%%%%%%%%%
% plot mean sacc by position (not time)
%%%%%%%%%
figure('name', 'mean saccade amplitude by position')
for b=1:5
    for a=5:-1:1
        clear S fixnums
        S= SaccsByPos{a,b};
        S(isnan(S(:,1)),:)= [];
        meansacc(a,b) = mean(S(:,1));
    end
end
bar3(meansacc)


%%%%%%%%%%%%
topleftSaccades = SaccsByPos{5,1}; topleftSaccades(isnan(topleftSaccades(:,1)),:)= [];
toprightSaccades = SaccsByPos{1,1}; toprightSaccades(isnan(toprightSaccades(:,1)),:)= [];
botleftSaccades = SaccsByPos{5,5};  botleftSaccades(isnan(botleftSaccades(:,1)),:)= [];
botrightSaccades = SaccsByPos{1,5};  botrightSaccades(isnan(botrightSaccades(:,1)),:)= [];
topSaccades = SaccsByPos{2:4,1};  topSaccades(isnan(topSaccades(:,1)),:)= [];
botSaccades =SaccsByPos{2:4,5};  botSaccades(isnan(botSaccades(:,1)),:)= [];
leftSaccades = SaccsByPos{5,2:4};  leftSaccades(isnan(leftSaccades(:,1)),:)= [];
rightSaccades = SaccsByPos{1,2:4};  rightSaccades(isnan(rightSaccades(:,1)),:)= [];
midSaccades = [SaccsByPos{2,2}; SaccsByPos{3,2}; SaccsByPos{4,2};...
    SaccsByPos{2,3}; SaccsByPos{3,3}; SaccsByPos{4,3};...
    SaccsByPos{2,4}; SaccsByPos{3,4}; SaccsByPos{4,4}];
midSaccades(isnan(midSaccades(:,1)),:)= [];


figure
subplot(3,3,1); rose(topleftSaccades(:,2));
subplot(3,3,2); rose(topSaccades(:,2));
subplot(3,3,3); rose(toprightSaccades(:,2));
subplot(3,3,4); rose(leftSaccades(:,2));
subplot(3,3,5); rose(midSaccades(:,2));
subplot(3,3,6); rose(rightSaccades(:,2));
subplot(3,3,7); rose(botleftSaccades(:,2));
subplot(3,3,8); rose(botSaccades(:,2));
subplot(3,3,9); rose(botrightSaccades(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flip corners around to get composite
figure('Position',[1 1 1000 1000])
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

figure('Position',[1 1 1200 1000])
% split up corner saccades by time
for t=1:3
    temp = [];
    subplot(4,4,t)
    for j=1:5
        temp = [temp; CornerTL(CornerTL(:,3)==(t-1)*5+j,1:2)]; %#ok<AGROW>
    end
    Corner(t).saccades = temp; %#ok<*AGROW>  
     csvwrite(['Corner_t=' int2str(t) '.csv'], Corner(t).saccades);
   map = SurfaceFromPoints(Corner(t).saccades);
      [Sx Sy] = pol2cart(temp(:,2), temp(:,1));      
    ContourPlotFromPoints([Sx Sy]);
   title(['Corner Saccades for t = ' int2str((t-1)*5+1) ' - ' int2str(t*5)]);
%    csvwrite(['Corner_t=' int2str(t) '.csv'],map);
end
temp = CornerTL(CornerTL(:,3)>15,1:2);
subplot(4,4,4)
Corner(4).saccades = temp;
 csvwrite(['Corner_t=' int2str(4) '.csv'], Corner(4).saccades);
map = SurfaceFromPoints(Corner(4).saccades);
[Sx Sy] = pol2cart(temp(:,2), temp(:,1));
ContourPlotFromPoints([Sx Sy]);
title(['Corner Saccades for t > ' int2str(15)]);
% csvwrite(['Corner_t=' int2str(4) '.csv'],map);
% split up EdgeL sacades by time

for t=1:3
    temp = [];
    subplot(4,4,t+4)
    for j=1:5
        temp = [temp; EdgeL(EdgeL(:,3)==(t-1)*5+j,1:2)]; %#ok<AGROW>
    end
    EdgeLeft(t).saccades = temp;
 csvwrite(['EdgeL_t=' int2str(t) '.csv'],EdgeLeft(t).saccades);
    map = SurfaceFromPoints(EdgeLeft(t).saccades);
       [Sx Sy] = pol2cart(temp(:,2), temp(:,1));      
    ContourPlotFromPoints([Sx Sy]);
%     csvwrite(['EdgeL_t=' int2str(t) '.csv'],map);
    title(['Vertical Edge Saccades for t = ' int2str((t-1)*5+1) ' - ' int2str(t*5)]);
end
temp = EdgeL(EdgeL(:,3)>15,1:2);
subplot(4,4,8)
EdgeLeft(4).saccades = temp;
 csvwrite(['EdgeL_t=' int2str(4) '.csv'],EdgeLeft(4).saccades);
map = SurfaceFromPoints(EdgeLeft(4).saccades);
   [Sx Sy] = pol2cart(temp(:,2), temp(:,1));      
    ContourPlotFromPoints([Sx Sy]);
title(['Vertical Edge Saccades for t > ' int2str(15)]);
% csvwrite(['EdgeL_t=' int2str(4) '.csv'],map);

% split up EdgeT saccades by time
for t=1:3
    temp = [];
    subplot(4,4,t+8)
    for j=1:5
        temp = [temp; EdgeT(EdgeT(:,3)==(t-1)*5+j,1:2)]; %#ok<AGROW>
    end
    EdgeTop(t).saccades = temp;
    csvwrite(['EdgeT_t=' int2str(t) '.csv'],EdgeTop(t).saccades);
    map = SurfaceFromPoints(EdgeTop(t).saccades);
       [Sx Sy] = pol2cart(temp(:,2), temp(:,1));      
    ContourPlotFromPoints([Sx Sy]);
    title(['Horizontal Edge Saccades for t = ' int2str((t-1)*5+1) ' - ' int2str(t*5)]);
%         csvwrite(['EdgeT_t=' int2str(t) '.csv'],map);
end
temp = EdgeT(EdgeT(:,3)>15,1:2);
subplot(4,4,4+8)
EdgeTop(4).saccades = temp; %#ok<*AGROW>
csvwrite(['EdgeT_t=' int2str(4) '.csv'],EdgeTop(4).saccades);
map = SurfaceFromPoints(EdgeTop(4).saccades);
   [Sx Sy] = pol2cart(temp(:,2), temp(:,1));      
    ContourPlotFromPoints([Sx Sy]);
title(['Horizontal Edge Saccades for t > ' int2str(15)]);
% csvwrite(['EdgeT_t=' int2str(4) '.csv'],map);

% and for mid saccades
for t=1:3
    temp = [];
    subplot(4,4,t+12)
    for j=1:5
        temp = [temp; midSaccades(midSaccades(:,3)==(t-1)*5+j,1:2)]; %#ok<AGROW>
    end
    middleSaccades(t).saccades = temp; %#ok<*AGROW>
   csvwrite(['midSaccades_t=' int2str(t) '.csv'],middleSaccades(t).saccades);
    map = SurfaceFromPoints(middleSaccades(t).saccades);
       [Sx Sy] = pol2cart(temp(:,2), temp(:,1));      
    ContourPlotFromPoints([Sx Sy]);
%     csvwrite(['midSaccades_t=' int2str(t) '.csv'],map);
    title(['Centre Region Saccades for t = ' int2str((t-1)*5+1) ' - ' int2str(t*5)]);
    
end
temp = midSaccades(midSaccades(:,3)>15,1:2);
subplot(4,4,4+12)
middleSaccades(4).saccades = temp; %#ok<*AGROW>
csvwrite(['midSaccades_t=' int2str(4) '.csv'],middleSaccades(4).saccades);
map= SurfaceFromPoints(middleSaccades(4).saccades);
   [Sx Sy] = pol2cart(temp(:,2), temp(:,1));      
    ContourPlotFromPoints([Sx Sy]);
% csvwrite(['midSaccades_t=' int2str(4) '.csv'],map);

title(['Centre Region Saccades for t > ' int2str(151)]);




    function S=fliphorizontally(S)
        S(:,2) = -(S(:,2)-pi/2)+pi/2;
        S(:,2) = mod(S(:,2),2*pi);
    end

    function S=flipvertically(S)
        S(:,2) = -(S(:,2)-pi/4)+pi/4;
        S(:,2) = S(:,2)+3*pi/2;
        S(:,2) = mod(S(:,2),2*pi);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below are my nested subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [a b] = getBin(x,y)
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
    end

    function Fixs = getFixesFromTextFile
        FixFile = fopen([Person int2str(run) 'FXD.txt']);
        data = fscanf(FixFile, '%c',inf);
        fclose(FixFile); clear FixFile
        %%% delete header
        cursor = regexp(data, 'Fix number', 'once');
        data(1:(cursor-1)) = [];
        cursor = regexp(data, '1', 'once');
        data(1:(cursor-1)) = [];
        % convert to matrix
        Fixs = str2num(data); %#ok<ST2NM>
    end

    function  Fixations = getFixInfo(Fixs)
        FinalFixMinDuration = 200;
        %define offset to turn Clearview output coords to screen coords.
        xoffset = (1200-1024)/2; yoffset = (1600-1024)/2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Note: I am including the final fixation on the fixation cross as
        % a trial data point... Otherwise I have no data for my first
        % saccade!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        index= find((Fixs(:,2)<=EndTime).*(Fixs(:,2)>=StartTime));
        if size(index,1)>0
            % sort out initial fixation, which probably started
            % before the onset of the stimulus. We will consider the time
            % from stimulus onset to end of fixation as a fixation. Unless
            % it lasts shorter than 200ms.
            if index(1)>1
                index = [index(1)-1; index];
            end
            CroppedInitFixLength = Fixs(index(1),2)+Fixs(index(1),3)-StartTime;
            if CroppedInitFixLength<200
                index(1) = [];
            end
            % now sort out final fixation. Only count if participant
            % fixates for at least 200ms before hitting button.
            TimeBetweenFinalFixAndButtonPress = EndTime- Fixs(index(size(index,1)),2);
            if  TimeBetweenFinalFixAndButtonPress < FinalFixMinDuration
                index(size(index,1)) = [];
            end
        end
        index= find((Fixs(:,2)<=EndTime).*(Fixs(:,2)>=StartTime));
        if size(index,1)==0
            % if no fixations recorded, assume time very short and there
            % was one fixation at hte fixation cross
            Fixations.number = 1;
        elseif size(index,1)>0
            Fixations.y = 1024-(1200-Fixs(index, 5)-xoffset);
            Fixations.x = Fixs(index, 4)-yoffset;
            Fixations.number= size(index,1);
        else
            % no fixations recorded for trial!!!
            Fixations.number = 0;
        end
    end % Fixations

    function [d theta] = GetSaccadeDistsAndAngles
        d = zeros(1, Fixations.number-1);
        theta = zeros(1, Fixations.number-1);
        for s  = 1:(Fixations.number-1)
            inBoundary = (Fixations.x(s)>0)&&(Fixations.x(s)<1025)...
                &&(Fixations.y(s)>0)&&(Fixations.y(s)<1025)...
                &&(Fixations.x(s+1)>0)&&(Fixations.x(s+1)<1025)...
                &&(Fixations.y(s+1)>0)&&(Fixations.y(s+1)<1025);
            if inBoundary==1
                xdiff = Fixations.x(s)-Fixations.x(s+1);
                ydiff = Fixations.y(s)-Fixations.y(s+1);
                d(s) = PixelsToVisualAngle(sqrt(xdiff^2+ydiff^2));
                if xdiff == 0
                    theta(s) = pi/2;
                elseif xdiff>0
                    theta(s) = atan(ydiff/xdiff);
                elseif xdiff <0
                    theta(s) = pi+atan(ydiff/xdiff);
                end
            else
                d(s)=NaN;
                theta(s)=NaN;
            end
        end
    end % SaccadeInfo
end

function Results = SortOutEventData(Person, numberofParts)
Results = zeros(360,13);
tick=0;
for erun = 1:numberofParts
    filename = [Person int2str(erun) 'EVD.txt'];
    EventFile = fopen(filename);
    data = fscanf(EventFile, '%c',inf);
    fclose(EventFile);
    %%% delete header
    cursor = regexp(data, 'Time', 'once');
    data(1:(cursor-1)) = [];
    N = size(data,2);
    % find fixcross.jpg.
    FixcrossLocations = regexp(data, 'fixcross.jpg');
    n = size(FixcrossLocations, 2);
    pngLocations = regexp(data, 'png');
    % for each fixcross
    % check time for following imageshow event (as it could be a rouge
    % keypress!)
    % then get start and stop time for trial.
    for t = 1:n
        cursor = FixcrossLocations(t);
        cursor = cursor + 13;
        followingevent = regexp(data(cursor:N), 'Keyboard|ShowSlide', 'once', 'match');
        % Catch any keypresses made while the Fixation Cross is
        % being displayed
        while strcmp(followingevent, 'Keyboard')
            cursor=regexp(data(cursor:N), 'Keyboard', 'once')+cursor;
            cursor=regexp(data(cursor+8:N), '[a-zA-Z]', 'once')+cursor;
            followingevent = regexp(data(cursor:N), 'Keyboard|ShowSlide', 'once', 'match');
        end
        if strcmp(followingevent, 'Keyboard')
            disp('bug to fix');
            break
        end
        StartTime = str2double(regexp(data(cursor:N), '[0-9][0-9][0-9][0-9]+', 'once', 'match'));
        Filename = regexp(data(cursor:N), 'RMS[._=\w]*png',  'once', 'match');
        cursor = pngLocations(t)+3;
        EndTime = str2double(regexp(data(cursor:N), '[0-9]*', 'once', 'match'));
        ReactionTime = (EndTime-StartTime)/1000;
        followingevent = regexp(data(cursor:N), 'Keyboard|ShowSlide', 'once', 'match');
        if strcmp(followingevent, 'ShowSlide')
            disp('bug to fix');
            break
        end
        Keypress =   str2double(regexp(data(cursor+10:N), '\s(13|32)\s', 'once', 'match') )  ;
        if Keypress==32
            Found = 1;
        else
            Found = 0;
            tick=tick+1;
        end
        RMS = str2double(regexp(Filename, '(?<=RMS=)(1(\.1)*)|0\.9(?=_)', 'match'));
        beta = str2double(regexp(Filename, '(?<=beta=)1\.[567]*', 'match'));
        seed = str2double(regexp(Filename, '(?<=seed=)[0-9]*', 'match'));
        r = str2double(regexp(Filename, '(?<=r=)[0-9]*', 'match'));
        Results(t+(erun-1)*90,1:8) = [StartTime, EndTime,  r, RMS, beta, seed,  ReactionTime, Found];
    end
end


end

function PData = creatEmptyResultsStructure
PData(1).Results = zeros(360,10);
PData(2).Results = zeros(360,10);
PData(3).Results = zeros(360,10);
PData(4).Results = zeros(360,10);
PData(5).Results = zeros(360,10);
PData(6).Results = zeros(360,10);
PData(7).Results = zeros(360,10);
end



