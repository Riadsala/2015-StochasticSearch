function [SaccadeDists SaccadeAngles] =AnalysisDataFromRoughnessExperiment
clear all close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 
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
People = {'mk'};%, 'my', 'lm', 'IHH', 'hw', 'am', 'ED'};
numberOfPeople = size(People, 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Empty Arrays to Save Results Into
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PData = creatEmptyResultsStructure;
NumberOfFixs = cell(3,3,2);
SaccadeDists = [];
SaccadeDistsP = cell(7);
SaccadeAngles = [];
IndvSaccadeAmps = cell(50,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Each Person, Extract Results from Clearview Export Text Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Pctr = 1:numberOfPeople
    Cov(Pctr).Vector = [];
    Person = People{Pctr};
    [PData(Pctr).Results timeoffsets] = SortOutEventData(Person, numberofParts);
    for run = 1:numberofParts
        run
        FixFile = fopen([Person int2str(run) 'FXD.txt']);
        data = fscanf(FixFile, '%c',inf);
        fclose(FixFile); clear FixFile
        %%% delete header
        cursor = regexp(data, 'Fix number', 'once');
        data(1:(cursor-1)) = [];
        cursor = regexp(data, '1', 'once');
        data(1:(cursor-1)) = [];
        % convert to matrix
        Fixs = str2num(data);
        for ctr = 1:90            
            t = ctr+(run-1)*90;
            % for each trial, get fix info
            StartTime = PData(Pctr).Results(t, 1)-timeoffsets(run);
            EndTime = PData(Pctr).Results(t,2)-timeoffsets(run);
            Fixations = getFixInfo;
            PData(Pctr).Results(t,9) = Fixations.number;
            [SaccAmps, SaccAngles] = GetSaccadeDistsAndAngles;
            if Fixations.number>1
                SaccadeDists = [SaccadeDists, SaccAmps];
                SaccadeAngles = [SaccadeAngles, SaccAngles];
                SaccadeDistsP{Pctr} = [SaccadeDistsP{Pctr}, SaccAmps];
            end
            for sa  = 1:50;                
                if Fixations.number>(sa+1)
                    CovVec(sa) = Fixations.c(sa);
                    IndvSaccadeAmps{sa} = [IndvSaccadeAmps{sa}, PixelsToVisualAngle(SaccAmps(sa))];
                else
                    CovVec(sa) = inf;
                end
            end
            Cov(Pctr).Vector= [ Cov(Pctr).Vector; CovVec];
            
        end
    end
    avgsac(Pctr) = size(SaccadeDistsP{Pctr},2)./(3*2*20*3);
    % get rough number of saccades taken before giving up
    m_index = PData(Pctr).Results(:,8)==0;
    % meanmiss(Pctr) = min((PData(Pctr).Results(m_index,7)))
    for beta  =1:3
        for RMS  =1:2
            for r = 1:3
                % Results(t,:) = [StartTime, EndTime,  r, RMS, beta, seed,...
                %   ReactionTime, FixNumber, FinalDist, medainDuration];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get index
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                index =  find((PData(Pctr).Results(:,3)==rSet(r))...
                    .*(PData(Pctr).Results(:,4)==RMSSet(RMS))...
                    .*(PData(Pctr).Results(:,5)==betaSet(beta))...
                    .*(PData(Pctr).Results(:,8)==1));
                % get number of fixation info
                 PData(Pctr).meanNumberofFixationsToTarget(beta, r, RMS) = mean(PData(Pctr).Results(index,9));
                meanNumberofFixationsToTarget(beta, r, RMS,Pctr) = mean(PData(Pctr).Results(index,9));
                PData(Pctr).StdDevNumberofFixationsToTarget(beta, r, RMS) = std(PData(Pctr).Results(index,9))./sqrt(size(index,1));
                NumberOfFixs{beta, r, RMS} = [NumberOfFixs{beta, r, RMS}; PData(Pctr).meanNumberofFixationsToTarget(beta, r, RMS)];
                %                 trialcoverage(beta,r,RMS, Pctr) = Fixations.c;
          Meancov(beta, r, RMS,Pctr) = mean(PData(Pctr).Results(index,10));
          MeancovRW(beta, r, RMS,Pctr) = mean(PData(Pctr).Results(index,11));
            end
        end
    end
    subplot(2,1,Pctr)
    plot(repmat(betaSet, [3,1])', simplesimmean, 'r')
end
save
subplotctr=0;

set(0,'DefaultAxesColorOrder',[0,0,0])
set(0,'DefaultAxesLineStyleOrder',{'-o','->','-s', '--o','-->','--s'})
% figure;
for beta  =1:3
    for r = 1:3
        subplotctr= subplotctr+1;
        subplot(3,3,subplotctr);
        mnf = meanNumberofFixationsToTarget(beta, r, RMS,:);
        
        mnf(1,1,1,8) = simplesimmean(beta,r);
        bar(mnf(:))
        title(['beta=' num2str(betaSet(beta)) ' and r=' int2str(rSet(r))]);
        %         outputforspss(:,subplotctr) = mnf(:);
    end
end
% csvwrite('outputforspss.csv', outputforspss);
meanNumFixs = zeros(3,3,1);
StdNumFixs = zeros(3,3,1);
for beta  =1:3
    for RMS  =1:2
        for r = 1:3
            meanNumFixs(beta,r,RMS) = mean(NumberOfFixs{beta, r, RMS});
            StdNumFixs(beta,r,RMS) = std(NumberOfFixs{beta, r, RMS});
        end
    end
end

meansa = zeros(1, 50);
stderrsa = zeros(1, 50);
 numtrials = zeros(1, 50);
for sa = 1:50
    meansa(sa) = mean(IndvSaccadeAmps{sa}) ;
    stderrsa(sa) = std(IndvSaccadeAmps{sa})./sqrt(size(IndvSaccadeAmps{sa},2)) ;
   numtrials(sa) = size(IndvSaccadeAmps{sa},2);
end
  wbltrends(1,:) = wblfit(IndvSaccadeAmps{1});
  
  for sa=2:50
%     wbltrends(sa,:) = wblfit([IndvSaccadeAmps{sa},IndvSaccadeAmps{sa+1},...
%         IndvSaccadeAmps{sa+2},IndvSaccadeAmps{sa+3},IndvSaccadeAmps{sa+4}]);
         wbltrends(sa,:) = wblfit(IndvSaccadeAmps{sa});
   
  end
figure
weibull_a = [5.3634, exp(1.8-0.0079*(1:49))];

weibull_b  =[1.7552, 1.8-0.0053*(1:49)];
subplot(2,2,1);
errorbar(meansa, stderrsa);
title('mean saccade amplitude');
xlabel('Saccade Number');
ylabel('Saccade Amplitude');
axis([0 55. 2.5 6]);
subplot(2,2,2);
plot(wbltrends(:,1),'b');
hold on
plot(weibull_a, 'r');
title('a paramter in weibull fit');
subplot(2,2,3);
plot(wbltrends(:,2),'b');
hold on
plot(weibull_b, 'r')
title('b paramter in weibull fit');
subplot(2,2,4);
plot(numtrials); 
title('number of trials');

% maxnumberofsaccs
% histogram saccade amplitudes

figure('Position',[1 1 800 400])
subplot(1,2,1)
SaccadeDists = PixelsToVisualAngle(SaccadeDists);
% index = find(SaccadeDists>20);
SaccadeDists(index) = [];
hist(SaccadeDists, 40, 'k');
% title('Saccade Amplitudes');
xlabel('Saccade Amplitude (\circ)', 'fontsize', 12);
ylabel('Number of Saccades', 'fontsize', 12)

% % rose plot saccade angles
subplot(1,2,2)
set(0,'DefaultAxesLineStyleOrder',{'-'})
rose(SaccadeAngles, 20);
% title('Saccade Directions')

set(0,'DefaultAxesLineStyleOrder',{'-o','->','-s'})


figure('Position',[1 1 800 400])
errorbar(repmat(betaSet, [3,1])',  meanNumFixs(:,:,2), StdNumFixs(:,:,2)/sqrt(numberOfPeople), 'b' );
hold all
%         title('RMS = 1.1', 'fontsize', 12);
%          errorbar(repmat(betaSet',1,3), ModelNumFixs(:, :, 2), ModelStdErr(:,:,2), 'r');
ylabel('Number of fixations', 'fontsize', 12);
axis([1.55, 1.75, 0, 45]);
plot(repmat(betaSet, [3,1])', simplesimmean, 'r')
xlabel('$\beta$', 'fontsize', 12, 'interpreter', 'latex'),


figure;
mc=  mean(Meancov,4);
sc = std(Meancov, [], 4);
errorbar([1.6, 1.65, 1.7], mc(:,1,2), sc(:,1,2)./sqrt(7), 'r')
hold on;
errorbar([1.6, 1.65, 1.7], mc(:,2,2), sc(:,2,2)./sqrt(7), 'b')
errorbar([1.6, 1.65, 1.7], mc(:,3,2), sc(:,3,2)./sqrt(7), 'c')
legend('r=100', 'r=225', 'r=350', 'location', 'southeast')
title('Human Trials');
xlabel('\beta');
ylabel('mean mean pixel-to-closest-fixation distance');

figure;
mc=  mean(MeancovRW,4);
sc = std(MeancovRW, [], 4);
errorbar([1.6, 1.65, 1.7], mc(:,1,2), sc(:,1,2)./sqrt(7), 'r')
hold on;
errorbar([1.6, 1.65, 1.7], mc(:,2,2), sc(:,2,2)./sqrt(7), 'b')
errorbar([1.6, 1.65, 1.7], mc(:,3,2), sc(:,3,2)./sqrt(7), 'c')
legend('r=100', 'r=225', 'r=350', 'location', 'southeast')
title('Random Walk');
xlabel('\beta');
ylabel('mean mean pixel-to-closest-fixation distance');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below are my nested subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    function  Fixations = getFixInfo
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
%             Fixations.cr = notSoSimpleSimf_justwalk(Fixations.number);
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
            d(s) = sqrt(xdiff^2+ydiff^2);
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


function [Results timeoffsets] = SortOutEventData(Person, numberofParts)
timeoffset = 0;
timeoffsets = [0 0 0];
Results = zeros(360,13);
tick=0;
for erun = 1:numberofParts
    timeoffsets(erun)= timeoffset;
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
        StartTime = StartTime+timeoffset;
        Filename = regexp(data(cursor:N), 'RMS[._=\w]*png',  'once', 'match');
        cursor = pngLocations(t)+3;
        EndTime = str2double(regexp(data(cursor:N), '[0-9]*', 'once', 'match'));
        EndTime = EndTime+timeoffset;
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
    timeoffset = timeoffset + EndTime;
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
