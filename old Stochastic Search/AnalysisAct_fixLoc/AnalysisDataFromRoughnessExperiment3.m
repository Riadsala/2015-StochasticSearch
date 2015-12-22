function AnalysisDataFromRoughnessExperiment2
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% .m file to read, process, and display results from Eye tracking
% experiment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params=getParams;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment Information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numberofParts = 4;
rSet = [100, 225, 350];
RMSSet = [0.9, 1.1];
betaSet = [1.6, 1.65, 1.7];
fixdistdist = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Participants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
People = {'lm', 'IHH', 'hw', 'am', 'ED', 'mk', 'my'};
numberOfPeople = size(People, 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Empty Arrays to Save Results Into
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PData = creatEmptyResultsStructure;

Hotspotmap = zeros(1024, 1024);
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
            beta = PData(Pctr).Results(trialnmbr,5);
            RMS = PData(Pctr).Results(trialnmbr,4);
            r=PData(Pctr).Results(trialnmbr,3);
            seed=PData(Pctr).Results(trialnmbr,6);
%             if (Fixations.number>1)
%                 filename = ['RMS=' num2str(RMS) '_beta=' num2str(beta) '_r=' int2str(r) '_seed=' int2str(seed) '.png']
%                 I = im2double(imread(filename));
%                 SalMap = AlasdairSaliencyMap(I, params);
%                 fixdist = AlasdairStochasticGazeCompFixLoc(I,SalMap, params, 1, 2, Fixations);
%                 fixdistdist = [fixdistdist fixdist];
%                 clear SalMap I 
%             end
        end % end trialctr loop
    end % end run loop   
end % end Pctr loop
 figure
 csvwrite('fixdist.csv', fixdistdist)
 hist(fixdistdist,50)
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

    function showHotSpotMap
       figure;
             
       disp(['number of fixations in hotspot map is ' int2str(sum(Hotspotmap(:)))]);
       Hotspotmapf = filter2(fspecial('gaussian',[99 99], 15), Hotspotmap, 'valid');
       imshow(Hotspotmapf,[]);
       size(Hotspotmapf)
 
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

function fixdist =  AlasdairStochasticGazeCompFixLoc(I,SalMap, params, seed1, seed2, Fixations)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alasdair's Saliency Alg
% Version Three - now checks itself if it has found the target
% Have also tidied up a lot.
%
% Version Two - Just do Gabor Filtering and max/median weighting
%
% What have I done?
% Removed contrast channel. The contrast channel is essentially a bank of
% bandpass filters. All the gabor channels for any given scale approximate
% a bandpass filter, so it loos like there is possible redundancy
% No point taking centre surrounds of orientation channel as that's jsut
% essentially applying a bandpass, and the gabor filters are already
% kind of like a bandpass filter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rand('state', seed1*seed2);
N=max(size(SalMap));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement Inhibition of Return (IOR) and Eccentricity Dependant
% Processing (EDP).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IOR_locations = zeros(params.IOR_decaytime,2); % this is used to store the locations where IOR is in effect.
IOR_mask = createIOR_mask;
% set initial fixation to image centre

TargetFixated = 0; numOfFixs = 0;
inputfilename = ['tmp'];
% now carry out search for target...
wrj = 0;
for t = 1:(Fixations.number-1)
    currentFixation = [Fixations.y(t), Fixations.x(t)];
    %                create inhibition of return map for this fixation
    IOR_map = createIOR_map;
    %                 apply exponential EDP mask and IOR mask
    tSalMap = SalMap.*EDPmask(currentFixation).*IOR_map;
    imwrite(SimplyNormalise(tSalMap), ['tSalMap' int2str(numOfFixs) '.png']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate new fixation - look at the 5 largest peaks in the saliency
    % map and choose one at random
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fSalMap = tSalMap;
    for f=1:3
        fixs(f,1:2) = max2(fSalMap);
        fixs(f,3) = tSalMap(fixs(f,1), fixs(f,2));
        fSalMap = fSalMap.*(1-circshift(IOR_mask, fixs(f,1:2)));
    end
%     imshow(tSalMap,[]);
%     hold all
%     plot([fixs(2,2),Fixations.x(1)], [fixs(2,1),Fixations.y(1)] ,'-bx')
%     plot([fixs(1,2),Fixations.x(1)], [fixs(1,1),Fixations.y(1)] ,'-bx')
%     plot([fixs(3,2),Fixations.x(1)], [fixs(3,1),Fixations.y(1)] ,'-bx')
%     plot(Fixations.x,Fixations.y, 'gx-')
    
    fixdist(t)=10000;
    
    for m=1:3
        fixdist(t) = min(fixdist(t), norm(fixs(m,1:2) - [Fixations.y(t+1), Fixations.x(t+1)]));
    end
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested Functions Live Below Here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function to create the IOR_mask
    function IOR_mask = createIOR_mask
        % use a gaussian distribution as a mask.
        IOR_mask = fspecial('gaussian',N, params.IOR_mask_size);
        % move gaussian to centre of image and normalise to [0, 1].
        IOR_mask = SimplyNormalise(circshift(IOR_mask, [-N/2 -N/2]));
    end
% function used to create the composite IOR map
    function IOR_map = createIOR_map
        % move all previous fixation locations back one step in time
        IOR_locations = circshift(IOR_locations, 1);
        % add current Fixation location to IOR list.
        IOR_locations(1, :) = currentFixation;
        %compute composite IOR map
        IOR_map = ones(N);
        for i = 1:min(numOfFixs, params.IOR_decaytime)
            IOR_map =IOR_map.*(1-(2/(i+1)).*circshift(IOR_mask, IOR_locations(i,:)));
        end
    end
% function used to create eccentricity dependant processing map.
    function EDP = EDPmask(Fixation)
        X = repmat((1:N)-Fixation(2), N,1);
        Y = repmat((1:N)'-Fixation(1),1,N);
        dist = sqrt(X.^2+Y.^2);
        EDP = exp(params.EDPradius.*dist);
    end
% function to check if the model is within params.FixationRadius of
% target's location
    function TargetFixated = CheckForTargetFixated(Fixation, Location)
        TargetFixatedx = (Fixation(1) > Location(1)-params.FixationRadius)&&...
            (Fixation(1) < Location(1)+params.FixationRadius);
        TargetFixatedy = (Fixation(2) > Location(2)-params.FixationRadius)&&...
            (Fixation(2) < Location(2)+params.FixationRadius);
        TargetFixated = TargetFixatedx * TargetFixatedy;
    end
end % end of main function


function DisplayOutputMap(map, loc, filename)
N = size(map,1);
map = SimplyNormalise(map);
loc1f = max(loc(1)-5,1);
loc1l = min(loc(1)+5,N);
loc2f = max(loc(2)-5,1);
loc2l = min(loc(2)+5,N);
output(:,:,1 ) = map;
output(:,:,2 ) = map;
output(:,:,3 ) = map;
output(loc1f:loc1l, loc2f:loc2l, 1)=output(loc1f:loc1l, loc2f:loc2l, 1)+0.5;
output(loc1f:loc1l, loc2f:loc2l, 2)=output(loc1f:loc1l, loc2f:loc2l, 1)-0.5;
output(loc1f:loc1l, loc2f:loc2l, 3)=output(loc1f:loc1l, loc2f:loc2l, 1)-0.5;
imwrite(output, filename);
end





