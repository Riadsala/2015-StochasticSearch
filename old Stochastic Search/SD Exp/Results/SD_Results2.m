function SD_Results2
close all
clear all
%define offset to turn Clearview output coords to image coords.
xoffset = (1200-1024)/2;
yoffset = (1600-1024)/2;
BlocksInRun = 27;
TrialsInBlock = 4;
TrialCounter = 0;
Results = zeros(2160, 5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sort out EventData - filenames and keypresses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for part = 1:7
    FixFile = fopen(['AC' int2str(part) 'FXD.txt']);
    FixData = fscanf(FixFile, '%c',inf);
    fclose(FixFile);
    %%% delete header
    cursor = regexp(FixData, 'Fix number', 'once');
    FixData(1:(cursor-1)) = [];
    cursor = regexp(FixData, '1', 'once');
    FixData(1:(cursor-1)) = [];
    % convert to matrix
    Fixs = str2num(FixData); %#ok<ST2NM>
    clear FixData
    EventFile = fopen(['AC' int2str(part) 'EVD.txt']);
    EventData = fscanf(EventFile, '%c',inf);
    fclose(EventFile);
    %%% delete header
    cursor = regexp(EventData, 'Description', 'once');
    EventData(1:(cursor+14))=[];
    N = size(EventData,2);
    % find fixcross.jpg.
    FixcrossLocations = regexp(EventData, 'fixcross.jpg');
    n = size(FixcrossLocations, 2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % go through all the trials in the run and get break times
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ImageLocations = regexp(EventData, '_beta');
    T = size(ImageLocations,2);
    BreakLocations = regexp(EventData, 'mask.jpg');
    BlockStartTime = str2double(regexp(EventData, '\d*', 'match', 'once'));
    for b = 1:BlocksInRun
        BlockStartTimes(b) = BlockStartTime;
        for t = 1:TrialsInBlock
            TrialCounter = TrialCounter+1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % find relvant fixation cross and get following image display time
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cursor = FixcrossLocations((b-1).*TrialsInBlock+t);
            TrialStartTime = regexp(EventData(cursor:N), '\d*', 'match', 'once');
            cursor = ImageLocations((b-1).*TrialsInBlock+t);    
            TrialEndTime = regexp(EventData(cursor:N), '\d\d\d\d*', 'match', 'once');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Check fixation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Fidx = find((Fixs(:,2)<TrialEndTime).*(Fixs(:,2)>TrialStartTimes(b)));
%         Fixations.y = 1024-(1200-Fixs(Fidx, 5)-xoffset);
%         Fixations.x = Fixs(Fidx, 4)-yoffset;
            
            cursor = ImageLocations((b-1).*TrialsInBlock+t);            
            % get surface parameters from filename
            Filename = regexp(EventData(cursor:N), 'beta[._=\w]*png',  'once', 'match');
            beta = str2double(regexp(Filename, '(?<=beta=)1\.[567]*', 'match'));
            r = str2double(regexp(Filename, '(?<=r=)(50|100|150|200|250|300|350|400)', 'match'));
            if r>0
                phi = str2double(regexp(Filename, '(?<=phi=)(0|45|90|135|180|225|270|315)', 'match'));
            else
                r = NaN;
                phi = NaN;
            end
            % get Right|Left Keyboard response
            Response =  regexp(EventData(cursor:N), '(Right)|(Left)',  'once', 'match');
            % check that response was found before the next trial was shown
            if (b-1)*TrialsInBlock+t~=T
                RespLoc =  regexp(EventData(cursor:N), '(Right)|(Left)',  'once');
                if (RespLoc+cursor)>ImageLocations((b-1).*TrialsInBlock+t+1)
                    disp('problem with missing response');
                    break
                end
            end
            Results(TrialCounter,1:3) = [beta, r, phi];
            if strcmp(Response, 'Right')
                Results(TrialCounter,4) = 0;
            elseif strcmp(Response, 'Left')
                Results(TrialCounter,4) = 1;
            end
        end
        cursor = cursor+regexp(EventData(cursor:N), 'Blank image',  'once');
        BlockEndTime = str2double(regexp(EventData(cursor:N), '\d*', 'match', 'once'));
        BlockStartTime = str2double(regexp(EventData(cursor:N), '\d*', 'match', 'once'));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check that fixation was mantained - centre of image
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Fidx = find((Fixs(:,2)<BlockEndTime).*(Fixs(:,2)>BlockStartTimes(b)));
        Fixations.y = 1024-(1200-Fixs(Fidx, 5)-xoffset);
        Fixations.x = Fixs(Fidx, 4)-yoffset;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sort out Results - how does roughness and r effect signal detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betaSet = [1.6, 1.65, 1.7];
rSet = 50:50:400;
AccuracyOverEcc = zeros(3,8);
TargetAbsentAccuracy = zeros(3,1);
TA2Plot = zeros(3,8);
for b = 1:3
    for r=1:8
        beta = betaSet(b);       
        index = find((Results(:,1)==beta).*(Results(:,2)==rSet(r)));
        AccuracyOverEcc(b,r) = mean(Results(index, 4));
    end
    index = find((Results(:,1)==beta).*isnan(Results(:,2)));
    TargetAbsentAccuracy(b) = 1-mean(Results(index, 4));
    TA2Plot(b,:) = TargetAbsentAccuracy(b)*ones(1,8);
end
plot([rSet; rSet; rSet]', 100.*AccuracyOverEcc', '-x')
legend('\beta=1.6', '\beta=1.65', '\beta=1.7')
hold on;
plot([rSet; rSet; rSet]',100.*TA2Plot', ':x')
xlabel('r (pixels)'); ylabel('% correct');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get overal TP accuracy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index = find(Results(:,2)>0);
mean(Results(index,4))