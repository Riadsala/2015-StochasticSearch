function Experiment2AFC

iLink.doILink = 1;


% set up experiment
addpath('commonFiles/')

subjNum = input('input subject number: ', 's');
dataFolder = strcat('Data/s', subjNum, '/');

dataFilename = strcat(dataFolder, 'results.mat');
iLink.edfdatafilename = strcat(dataFolder, 'results.edf');

% define noise levels for easy and hard
%A = getContrastLevels(subj);
beta = [1.6, 1.65, 1.7];
N = 1024;

try
    load(dataFilename);   
catch
    nAngles = 8;
    nEccs = 4;
    nBlocks = 2;
    mkdir(dataFolder);
    % make conditions
    ctr = 0;
    for blk = 1:nBlocks
        for b = 1:length(beta);
            for phi = 1:nAngles
                for ecc = 1:nEccs
                    ctr = ctr + 1;
                    dat(ctr).beta = beta(b);
                    dat(ctr).phi = (phi/nAngles) * (2*pi);
                    dat(ctr).ecc = ecc*  N/(2*(nEccs+1));
                    dat(ctr).res = [];
                end
            end
            ctr = ctr  +1;
            dat(ctr).beta = beta(b);
            dat(ctr).phi = 0;
            dat(ctr).ecc = 0;
            dat(ctr).res = [];
        end
    end
end

blocksPerSession = 6;
fixCrossDisplayTime = 1.00;
stimuliDisplayTime = 0.25;
interStimInterval = 0.75;
penaltyDisplayTime = 2.00;

numTrials = 20;

%% set things up
bkgrndGreyLevel = 127;

[stimuliScrn] = Screen('OpenWindow',1, bkgrndGreyLevel);%, [1001 301 1600 900]
[display.width, display.height]=Screen('WindowSize', stimuliScrn);%screen returns the size of the window
display.midX = round(display.width/2);
display.midY = round(display.height/2);
% HideCursor;
% blank screen
blank = bkgrndGreyLevel*ones(N);
t_blank = Screen('MakeTexture', stimuliScrn, blank);
blank(round(N/2-3):round(N/2+3),round(N/2-3):round(N/2+3)) = 64;
t_fixDot = Screen('MakeTexture', stimuliScrn, blank);
clear blank
% penalty screen
red = 255*ones(N,N,3);
red(:,:,2:3) = 0;
t_red = Screen('MakeTexture', stimuliScrn, red);
clear red;
% fixation cross
fixCross = makeFixCross(N, bkgrndGreyLevel);
t_fixCross = Screen('MakeTexture', stimuliScrn, fixCross);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize all eyelink and calibrate
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iLink.doILink == 1
    iLink = InitEyeLink(iLink, stimuliScrn);
    iLink.howLongToTryFor = 3;
    iLink.centralFixDur = fixCrossDisplayTime;
    iLink.stimuliDisplayTime = stimuliDisplayTime;
end

for BLOCK = 1:blocksPerSession
    % randomly select a block to work on
    blk  = pickABlock2Do;
    % create results array
    dat(blk).res = zeros(numTrials,1);
    
    % if we're not doing the inital block, recalibrate
    if (BLOCK > 1) && (iLink.doILink == 1)
        EyelinkDoTrackerSetup(iLink.el);
    end
    WaitSecs(0.25);
    % draw start block screen screen
    Screen('DrawTexture',stimuliScrn,t_blank);
    message = (['About to start block. \n You will first see preview of target location. \n When happy, press a key']); ;
    DrawFormattedText(stimuliScrn, message, 'center', 'center', 200);
    Screen('Flip',stimuliScrn,0);
    getObserverInputFS;
    
    % draw target preview
    [x,y] = pol2cart(dat(blk).phi, dat(blk).ecc);
                 x = round(x+N/2);
                y = round(y+N/2);
    
    
    stimPreview  = bkgrndGreyLevel*ones(N);
    stimPreview((x-10):(x+10), (y-10):(y+10)) = 200;
    t_preview =  Screen('MakeTexture',stimuliScrn,stimPreview);
    Screen('DrawTexture',stimuliScrn,t_preview);
    Screen('Flip',stimuliScrn,0);
    WaitSecs(0.5);
    getObserverInputFS;
    
    % do trials
    t = 0;
    while t < numTrials
        resp = carryOutATrial(dat(blk), t<numTrials);
        
        if resp == 3
            sca;
            ShowCursor;
            return;
        elseif resp == -1
            % failed to get stable central fixation - recalibrate
            EyelinkDoTrackerSetup(iLink.el);
            % now remind observer where target will appear
            % draw target preview
            stimPreview((x-10):(x+10), (y-10):(y+10)) = 200;
             t_preview =  Screen('MakeTexture',stimuliScrn,stimPreview);
             Screen('DrawTexture',stimuliScrn,t_preview);
            Screen('Flip',stimuliScrn,0);
            WaitSecs(0.5);
            getObserverInputFS;
            
        elseif resp == -2
            % user broke central fixation - discard trial
            % (do nothing)
        else
            t = t + 1;
            dat(blk).res(t) = resp;
        end
        
    end
    Screen('Close', t_preview)
    save(dataFilename, 'dat');
end

 Screen('Close', stimuliScrn);
ShowCursor;

save(dataFilename, 'dat');

%% Nested Functions

    function resp = carryOutATrial(dat,showPostPreview)
        
        % randomly select a TA trial  
        z = randperm(10);
        stimTA  = imread(strcat('ims/ta_beta', num2str(dat.beta), '_', int2str(z(1)), '.png'));
        % randomly select a TP trial
        stimTP  = imread(strcat('ims/tp_beta=', num2str(dat.beta), ...
        '_r=', int2str(dat.ecc), '_phi=', num2str(180*dat.phi/pi), '_', int2str(z(2)), '.png'));
        
        %% a trial
        if rand<0.5;
            targetint= 0;
            stim{1} = stimTP;
            stim{2}= stimTA;
        else
            targetint = 1;
            stim{1} = stimTA;
            stim{2}= stimTP;
        end
        t_stimulus1 = Screen('MakeTexture',stimuliScrn,stim{1});
        t_stimulus2 = Screen('MakeTexture',stimuliScrn,stim{2});
        
        %% wait for central fixation before showing stimuli
        Screen('DrawTexture', stimuliScrn, t_fixCross);
        Screen('Flip',stimuliScrn,0);
        
        if iLink.doILink
            % Start data recording to EDF file, BEFORE DISPLAY. */
            Eyelink('startrecording');
            % record a few samples before we actually start displaying
            WaitSecs(0.1);
            Eyelink('message', 'DISPLAY_ON');	 % message for RT recording in analysis
            Eyelink('message','SYNCTIME');
            % WAIT for central fixation
            status = WaitForCentralFixation(iLink, display);
            if status <0
                resp = -1;
                return
            end
        else
            % not using eyelink, so do simple version
            WaitSecs(fixCrossDisplayTime);
            status = 1;
        end
        
        %% display stimulus
        Screen('DrawTexture',stimuliScrn,t_stimulus1);
        Screen('Flip',stimuliScrn,0);
        if iLink.doILink
            keptCentralFixation = KeepCentralFixation(iLink, display);
        else
            WaitSecs(stimuliDisplayTime)
            keptCentralFixation = true;
        end
        % only continue with trial if central fixation was matained
        if keptCentralFixation
            % display fixation dot
            Screen('DrawTexture',stimuliScrn,t_fixCross);
            Screen('Flip',stimuliScrn,0);
            if iLink.doILink
                keptCentralFixation = KeepCentralFixation(iLink, display, interStimInterval);
            else
                WaitSecs(interStimInterval)
                keptCentralFixation = true;
            end
        end
        % only continue with trial if central fixation was matained
        if keptCentralFixation
            % display stimulus
            Screen('DrawTexture',stimuliScrn,t_stimulus2);
            Screen('Flip',stimuliScrn,0);
            if iLink.doILink
                keptCentralFixation = KeepCentralFixation(iLink, display);
            else
                WaitSecs(stimuliDisplayTime)
                keptCentralFixation = true;
            end
        end
        
        %% get response/punish
        if keptCentralFixation
            %% trial good, so display blank and wait for response.
            Screen('DrawTexture',stimuliScrn,t_blank);
            Screen('Flip',stimuliScrn,0);
            % wait for response key
            responseKeyHit = 0;
            while ~responseKeyHit
                [resp responseKeyHit] = getObserverInputFS;
                WaitSecs(0.001);
                resp = resp==targetint;
            end
            
        else
            %% observer broke central fixation, so show anrgy screen
            Screen('DrawTexture',stimuliScrn,t_red);
            Screen('Flip',stimuliScrn,0);
            WaitSecs(penaltyDisplayTime);
            resp = -2;
        end
        
        % show blank
        Screen('DrawTexture',stimuliScrn,t_blank);
        Screen('Flip',stimuliScrn,0);
        WaitSecs(0.1);
        
        % post response cue
        if showPostPreview
            cue = bkgrndGreyLevel * ones(N);
            cue((x-5):(x+5), (y-5):(y+5)) = 65;
            t_stimulus = Screen('MakeTexture',stimuliScrn,cue);
            Screen('DrawTexture',stimuliScrn,t_stimulus);
            Screen('Flip',stimuliScrn,0);
            WaitSecs(0.5);
            Screen('Close', t_stimulus)
        end
        Screen('Close', t_stimulus1)
        Screen('Close', t_stimulus2)
        if iLink.doILink
            Eyelink('stoprecording');
        end
        
    end

    function blk  = pickABlock2Do
        todo = [];
        for c = 1:length(dat)
            if isempty(dat(c).res)
                todo = [todo, c]; %#ok<AGROW>
            end
        end
        disp(strcat('Number of blocks left: ', int2str(length(todo))));
        blk = todo(randi(length(todo)));
    end

end


