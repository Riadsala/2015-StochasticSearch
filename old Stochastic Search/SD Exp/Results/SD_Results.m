function SD_Results
figure
%define offset to turn Clearview output coords to image coords.
xoffset = (1200-1024)/2;
yoffset = (1600-1024)/2;

BlocksInRun = 16;
TrialsInBlock = 18;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sort out EventData - filenames and keypresses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for part = 1:2
    EventFile = fopen(['AC' int2str(part) 'EVD.txt']);
    EventData = fscanf(EventFile, '%c',inf);
    fclose(EventFile);
    %%% delete header
    cursor = regexp(EventData, 'Time', 'once');
    EventData(1:(cursor-1)) = [];
    N = size(EventData,2);
    % find fixcross.jpg.
    FixcrossLocations = regexp(EventData, 'fixcross.jpg');
    n = size(FixcrossLocations, 2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % go through all the trials in the run and get break times
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BreakLocations = regexp(EventData, 'mask.jpg');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % go through all the trials in the run and get keypresses
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ImageLocations = regexp(EventData, '_beta');
    T = size(ImageLocations,2);
    for t = 1:T
        pt = (part-1).*288+t;
        % move cursor to correct location
        cursor = ImageLocations(t);
        % get surface parameters from filename
        Filename = regexp(EventData(cursor:N), 'beta[._=\w]*png',  'once', 'match');
        beta = str2double(regexp(Filename, '(?<=beta=)1\.[567]*', 'match'));
        r = str2double(regexp(Filename, '(?<=r=)(60|120|180|240|300|360|420)', 'match'));
        if r>0
            phi = str2double(regexp(Filename, '(?<=phi=)(0|45|90|135|180|225|270|315)', 'match'));
        else
            r = NaN;
            phi = NaN;
        end
        % get Right|Left Keyboard response
        Response =  regexp(EventData(cursor:N), '(Right)|(Left)',  'once', 'match');
        % check that response was found before the next trial was shown
        if t~=T
            RespLoc =  regexp(EventData(cursor:N), '(Right)|(Left)',  'once');
            if (RespLoc+cursor)>ImageLocations(t+1)
                disp('problem with missing response');
                t
                break
            end
        end
        Results(pt,1:3) = [beta, r, phi];
        if strcmp(Response, 'Right')
            Results(pt,4) = 0;
        elseif strcmp(Response, 'Left')
            Results(pt,4) = 1;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % look at fixations and determine whether participant mvoed from
    % fixation cross when they shouldn't of
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FixFile = fopen(['AC' int2str(part) 'FXD.txt']);
    FixData = fscanf(FixFile, '%c',inf);
    fclose(FixFile);
    %%% delete header
    cursor = regexp(FixData, 'Fix number', 'once');
    FixData(1:(cursor-1)) = [];
    cursor = regexp(FixData, '1', 'once');
    FixData(1:(cursor-1)) = [];
    % convert to matrix
    Fixs = str2num(FixData);
    clear FixData
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sort out Results - how does roughness and r effect signal detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betaSet = [1.6, 1.65, 1.7];
rSet = [60, 120, 180, 240, 300, 360, 420];
AccuracyOverEcc = zeros(3,6);
TargetAbsentAccuracy = zeros(3,1);
TA2Plot = zeros(3,7);
for b = 1:3
    for r=2:7
        if b==1
            beta = betaSet(b);
            index = find((Results(:,1)==beta).*(Results(:,2)==rSet(r-1)));
            AccuracyOverEcc(b,r-1) = mean(Results(index, 4));
            AccuracyOverEcc(b,7) = NaN;
        else
            beta = betaSet(b);
            index = find((Results(:,1)==beta).*(Results(:,2)==rSet(r)));
            AccuracyOverEcc(b,r) = mean(Results(index, 4));
            AccuracyOverEcc(b,1) = NaN;
        end
    end
    index = find((Results(:,1)==beta).*isnan(Results(:,2)));
    TargetAbsentAccuracy(b) = 1-mean(Results(index, 4));
    TA2Plot(b,:) = TargetAbsentAccuracy(b)*ones(1,7);
end
plot([rSet; rSet; rSet]', 100.*AccuracyOverEcc', '-x')
legend('\beta=1.6', '\beta=1.65', '\beta=1.7')
hold on;
plot([rSet; rSet; rSet]',100.*TA2Plot', ':x')
xlabel('r (pixels)'); ylabel('% correct');