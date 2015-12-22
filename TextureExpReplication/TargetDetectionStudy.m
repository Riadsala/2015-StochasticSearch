function VisualSearchStudy
addpath('commonfiles/');
participantNumber = input('enter participant number:');

iLink.doILink = 1;

%% set things up
bkgrndGreyLevel = 127;
N = 512;
[stimuliScrn] = Screen('OpenWindow',1, bkgrndGreyLevel);%, [001 01 1600 900]
[params.width, params.height]=Screen('WindowSize', stimuliScrn);%screen returns the size of the window
params.midX = round(params.width/2);
params.midY = round(params.height/2);

HideCursor;

params.fixationCrossDuration = 1;
params.maxStimuliparams = 0.200;

params.penaltyTime = 3;
params.postTrialTime = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize all eyelink and calibrate
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iLink.doILink == 1
    iLink.edfdatafilename = strcat('ac_', int2str(participantNumber), '_td');
    iLink = InitEyeLink(iLink, stimuliScrn);
    iLink.howLongToTryFor = 5;
    iLink.centralFixDur = params.fixationCrossDuration;
    iLink.centralTol = 100^2;
end

%get filenames
trials= dir('im/detection/*.png');
for t= 1:length(trials)
    trials(t).beta =  regexp(trials(t).name, '(?<=beta=)\d\.[\d]*', 'match');
    isTP = isempty(regexp(trials(t).name, 'TA', 'match'));
    trials(t).seed =  regexp(trials(t).name, '(?<=seed=)[\d]*', 'match');
    if isTP
        trials(t).isTP = 1;
        trials(t).r =  regexp(trials(t).name, '(?<=r=)[\d]*', 'match');
        trials(t).phi =  regexp(trials(t).name, '(?<=phi=)[\d]*', 'match');
        
    else
        trials(t).r = nan;
        trials(t).phi = nan;
        trials(t).isTP = 0;
    end
end

% randomly permute
trials = trials(randperm(length(trials)));
practrials = trials(1:20);
trials = trials(randperm(length(trials)));
%% params intro screen
setupBlankAndFixCross;
Screen('DrawTexture', stimuliScrn, params.t_blank);
DrawFormattedText(stimuliScrn, '20 practise trials - Press space when ready to start', 'center', 'center')
Screen('Flip', stimuliScrn)
WaitSecs(0.5);
KbWait;
GetSecs; % run to load function into memory

%% practise trials
for t= 1:length(practrials)
    resp = DoATrial(trials(t), stimuliScrn, params, iLink);    
end

Screen('DrawTexture', stimuliScrn, params.t_blank);
DrawFormattedText(stimuliScrn, 'Press space when ready to start', 'center', 'center')
Screen('Flip', stimuliScrn)
WaitSecs(0.5);
KbWait;

fout = fopen(strcat('results/ac_', int2str(participantNumber), '_td.txt'), 'w');
fprintf(fout, 'participantNumber, trialNum, TargetPresent, beta, seed, r, phi, respose\n');
for t= 1:length(trials)
    if mod(t,50)==0
        % time for a break
        Screen('DrawTexture', stimuliScrn, params.t_blank);
        DrawFormattedText(stimuliScrn, ['Time for a break\n ' int2str(t) ' out of ' int2str(length(trials))], 'center', 'center')
        Screen('Flip', stimuliScrn)
        WaitSecs(5);
        KbWait;
        EyelinkDoTrackerSetup(iLink.el);
    end
    resp = DoATrial(trials(t), stimuliScrn, params, iLink);
    trials(t)
    if trials(t).isTP
        fprintf(fout, '%d, %d, %d, %s, %s, %s, %s, %d\n', participantNumber, t,...
            trials(t).isTP, trials(t).beta{1}, trials(t).seed{1}, trials(t).r{1}, trials(t).phi{1}, resp);
    else
        fprintf(fout, '%d, %d, %d, %s, %s, NaN, NaN, %d\n', participantNumber, t,...
            trials(t).isTP, trials(t).beta{1}, trials(t).seed{1}, resp);
    end
    
end
fclose(fout);
Screen('DrawTexture', stimuliScrn, params.t_blank);
DrawFormattedText(stimuliScrn, 'Thank You!', 'center', 'center')
Screen('Flip', stimuliScrn)


if iLink.doILink
    EYELINK('closefile');
    EYELINK('shutdown');
end


KbWait;

sca
ShowCursor;
end

function resp = DoATrial(trial, stimuliScrn, params, iLink)
status=Eyelink('startrecording');
% load image
im = imread(['im/detection/' trial.name]);
t_im = Screen('MakeTexture', stimuliScrn, im);


% show fixation cross
Screen('DrawTexture', stimuliScrn, params.t_fixCross);
Screen('Flip', stimuliScrn)
% WaitSecs(params.fixationCrossDuration);
out = -1;
while out == -1
    out = WaitForCentralFixation(iLink, params);
    if out == -1
        Eyelink('stoprecording');
        EyelinkDoTrackerSetup(iLink.el);
        status=Eyelink('startrecording');
        % show fixation cross
        Screen('DrawTexture', stimuliScrn, params.t_fixCross);
        Screen('Flip', stimuliScrn)
    end
end



% show stimulus
Screen('DrawTexture', stimuliScrn, t_im);
Screen('Flip', stimuliScrn);

status = KeepCentralFixation(iLink, params, params.maxStimuliparams);
if status
    Screen('DrawTexture', stimuliScrn, params.t_mask);
    Screen('Flip', stimuliScrn);
    WaitSecs(0.5);
    Screen('DrawTexture', stimuliScrn, params.t_blank);
    Screen('Flip', stimuliScrn);
    resp = getObserverInputYesNo;
else
    Screen('DrawTexture', stimuliScrn, params.t_red);
    Screen('Flip', stimuliScrn);
    WaitSecs(2);
    resp = -1;
end
Eyelink('stoprecording');
Screen('Close', t_im);
end
