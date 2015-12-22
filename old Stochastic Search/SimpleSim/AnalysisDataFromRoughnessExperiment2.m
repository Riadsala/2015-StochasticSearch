function AnalysisDataFromRoughnessExperiment2
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ar
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
IndvSaccadeAmps20 = cell(20,1);
CoverageTimeSeries = zeros(50, 360, 7);
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
            SaccadeAngles = [SaccadeAngles, SaccAngles];
            SaccadeDists = [SaccadeDists, SaccDists];
%             CoverageVector = zeros(50,1);
            for f = 1:49;
                if Fixations.number>=(f+2)
%                     CoverageVector(f) = Fixations.c(f);
                    IndvSaccadeAmps{f} = [IndvSaccadeAmps{f}, SaccDists(f)];
%                    else
%                     CoverageVector(f) = inf;
                end
            end
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % Get first 20 saccades from all trials which had more than 25
%             % fixations
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if Fixations.number>24
%                 for f = 1:20
%                     IndvSaccadeAmps20{f} = [IndvSaccadeAmps20{f}, SaccDists(f)];
%                 end
%             end
%             CoverageTimeSeries(:,trialnmbr, Pctr) = CoverageVector;
        end % end trialctr loop
    end % end run loop
end % end Pctr loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get basic stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interPerson_meanNumFix = zeros(3,3,2);
interPerson_sderNumFix = zeros(3,3,2);
for beta  =1:3
    for RMS  =1:2
        for r = 1:3
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get index of relevant files
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            index =  find(...
                (PData(Pctr).Results(:,3)==rSet(r))...
                .*(PData(Pctr).Results(:,4)==RMSSet(RMS))...
                .*(PData(Pctr).Results(:,5)==betaSet(beta)));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get number of fixations to target for each person
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get number of fixation info
            NumberofFixationsToTarget = zeros(7,1);
            for Pctr = 1:numberOfPeople
                NumberofFixationsToTarget(Pctr) = mean(PData(Pctr).Results(index,9));                
            end
            interPerson_meanNumFix(beta, r, RMS) = mean(NumberofFixationsToTarget);
            interPerson_sderNumFix(beta, r, RMS) = std(NumberofFixationsToTarget)./sqrt(numberOfPeople);
        indiv_meanNumFix(beta,r,RMS,:) = NumberofFixationsToTarget;
        end
    end
end



save

plotNumFix % plot effect of beta and r on number of fixations to target
plotSacc_t % plot time series for saccade amplitudes and Weibull fits
% plotSaccHist
plotSaccHistAndRose % show overall saccade amplitudes and directions
% plotCoverage
plotModelHumanComp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below are my nested subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    function plotModelHumanComp
        model = csvread('meansacc.csv');
        figure
        plotctr = 0;
       for beta=1:3
           for r=1:3
               plotctr= plotctr + 1;
              subplot(3,3,plotctr);
              data(1:7) = indiv_meanNumFix(beta,r,2,:);
              data(8) = model(beta,r);
              bar(data)
              title(['\beta = ' num2str(betaSet(beta)) ' and r = ' int2str(rSet(r))]) ;
           end
       end
    end

    function plotNumFix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot graphs - number of fixations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        beta = repmat(betaSet, [3,1])';
        figure;
        errorbar(beta, interPerson_meanNumFix(:,:,2), interPerson_sderNumFix(:,:,2));
        xlabel('\beta');
        ylabel('interpersonal mean number of fixations to target');
        title('\beta, r and number of fixations to target (\sigma_(RMS)=1.1 only');
        hold on
         model = csvread('meansacc.csv');
         plot(model)
        
    end

    function plotSacc_t
        figure
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot graphs - saccade amplitude time-series
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        meanSaccadeAmp = zeros(1, 49);
        sderSaccadeAmp = zeros(1, 49);
        numtrialsfor_t = zeros(1, 49);
        weibullltrends = zeros(49, 2);
        for t = 1:49
            meanSaccadeAmp(t) = mean(IndvSaccadeAmps{t}) ;
            sderSaccadeAmp(t) = std(IndvSaccadeAmps{t})./sqrt(size(IndvSaccadeAmps{t},2)) ;
            numtrialsfor_t(t) = size(IndvSaccadeAmps{t},2);
%             weibullltrends(t,:) = wblfit(IndvSaccadeAmps{t});
        end
        % compute best fit regression lines for Weibull paramters
        weibull_a = [5.3634, exp(1.8-0.0079*(1:49))];
        weibull_b  =[1.7552, 1.8-0.0053*(1:49)];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot Saccade Amplitude Time-Series
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         subplot(2,1,1);
errorbar(meanSaccadeAmp, sderSaccadeAmp);
        title('Effect of Time on mean Saccade Amplitude');
        xlabel('Saccade Number'); ylabel('Saccade Amplitude');
        axis([0 55. 2.5 6]);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % plot number of trials in each data point
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         subplot(2,2,2); plot(numtrialsfor_t);
%         title('number of trials in each data point');
%         xlabel('Saccade Number'); ylabel('Number of trials in data points');
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % plot how time changes Weibull parameters
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         subplot(2,2,3);
%         plot(weibullltrends(:,1),'b');
%         hold on;
%         plot(weibull_a, 'r');
%         title('a paramter in weibull fit');
%         subplot(2,2,4);
%         plot(weibullltrends(:,2),'b');
%         hold on
%         plot(weibull_b, 'r')
%         title('b paramter in weibull fit');
    end

    function plotSaccHistAndRose
        figure('Position',[1 1 800 400])
        subplot(1,2,1)
        hist(SaccadeDists, 0.25:0.5:24.75, 'k');
        title('Saccade Amplitudes');
        xlabel('Saccade Amplitude (\circ)', 'fontsize', 12);
        ylabel('Number of Saccades', 'fontsize', 12)
        subplot(1,2,2)
        set(0,'DefaultAxesLineStyleOrder',{'-'})
        rose(SaccadeAngles, 20);
        title('Saccade Directions')
    end

    function plotCoverage
        figure; hold on;
        meanCoverage  =zeros(7,50);
        for p=1:7            
            for t=1:49
                idx = isfinite(CoverageTimeSeries(t,:,p));
                meanCoverage(p,t) = mean(CoverageTimeSeries(t,idx,p));
            end
            
        end
        plot(log(mean(meanCoverage,1)), 'b');
        axis([0 50, -2 5]);
        xlabel('Fixation Number');
        ylabel('log Coverage reduced by...');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot random walk results
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hold on
        trialcoverage = csvread('randomwalktrialcoverage.csv');
        plot(1:49, log(mean(trialcoverage)), 'k:');
         trialcoverage = csvread('randomwalk2trialcoverage.csv');
        plot(1:49, log(mean(trialcoverage)), 'g-.');
        a = csvread('randomcoordstrialcoverage.csv');
        plot(1:49, log(mean(a)), 'r-');
        legend('mean human performance', 'modelled based random walk', 'empirically based random walk', 'random coordinates');
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
            Fixations.c = coverage2([512 512], 1024);            
        elseif size(index,1)>0
            Fixations.y = 1024-(1200-Fixs(index, 5)-xoffset);
            Fixations.x = Fixs(index, 4)-yoffset;
            Fixations.number= size(index,1);
            Fixations.c = coverage2([Fixations.x, Fixations.y], 1024);            
        else
            % no fixations recorded for trial!!!
            Fixations.number = 0;
            Fixations.c = inf;
            Fixations.cr = inf;
            
        end
    end % Fixations

    function [d theta] = GetSaccadeDistsAndAngles
        d = zeros(1, Fixations.number-1);
        theta = zeros(1, Fixations.number-1);
        
        for s  = 1:(Fixations.number-1)
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
