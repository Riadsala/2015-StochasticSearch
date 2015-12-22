function EstimateSaccDist
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
numberOfPeople = length(People);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Empty Arrays to Save Results Into
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PData = creatEmptyResultsStructure;


Q = 32;
SaccByPos = zeros(Q, Q, Q, Q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Each Person, Extract Results from Clearview Export Text Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fout = fopen('clarke2009Fixations.txt', 'w');
fprintf(fout, 'Obs t f x y\n')
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
             if Fixations.number >1
            for f = 1:Fixations.number
               fprintf(fout, '%f %f %f %f %f\n', Pctr, trialnmbr, f, Fixations.x(f), Fixations.y(f));
             end
            end
            PData(Pctr).Results(trialnmbr,9) = Fixations.number;
            sacc = GetSaccadeDistsAndAngles;

            % go through each saccade and bin it by fixation location
            sacc = ceil(sacc./Q);

            % don't use last saccade, as that's probably towards the target
            for s = 1:(size(sacc,1)-1)
               SaccByPos(sacc(s,1),sacc(s,2),sacc(s,3),sacc(s,4)) = ...
                   SaccByPos(sacc(s,1),sacc(s,2),sacc(s,3),sacc(s,4)) + 1;
                
            end
        end % end trialctr loop
    end % end run loop
end % end Pctr loop
fclose(fout);
%% now want to do some smoothing
% create 4D filter!
F = smoother4DFilter(11, 3);
SaccByPos = convn(SaccByPos, F, 'same');

save SaccDistribution SaccByPos
% %% now we fold! [just the fixation locations, not the saccade targets though
% SaccByPos(:,1:(Q),:,:) = SaccByPos(:,1:(Q),:,:)+SaccByPos(:,(Q):-1:1,:,:);
% SaccByPos(1:(Q),:,:,:) = SaccByPos(1:(Q),:,:,:)+SaccByPos((Q):-1:1,:,:,:);
% 
% slice(:,:) = SaccByPos(14,14,:,:); imshow(slice,[])

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

    function [toFrom] = GetSaccadeDistsAndAngles
        toFrom = zeros(Fixations.number-1, 4);
        for s  = 1:(Fixations.number-1)
            toFrom(s,:) = [Fixations.x(s), Fixations.y(s), ...
                Fixations.x(s+1), Fixations.y(s+1)];
        end
        
        % remove saccades that start or end outside of stimuli
        toFrom(sum(toFrom>1024,2)>0,:) = [];
        toFrom(sum(toFrom<1,    2)>0,:) = [];
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


function output = PixelsToVisualAngle(input)

% input is a distance given in pixels

pixelsize = 0.255; % in mm
distfromdisplay = 0.87; % in meters

saccadelength = (pixelsize/1000)*input;

output = atand(saccadelength./distfromdisplay);

end
