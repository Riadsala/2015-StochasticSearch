function VisualSearchStudy2
addpath('commonfiles/');
participantNumber = input('enter participant number:');
personDat.id = participantNumber;
personDat.age = input('enter age:');
personDat.gender = input('enter gender:', 's');
save(strcat('results/personDat', int2str(participantNumber), '.mat'), 'personDat')

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
params.maxStimuliparams = 60;
params.catchTrialparams = [30, 15, 5];
params.penaltyTime = 3;
params.postTrialTime = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize all eyelink and calibrate
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iLink.doILink == 1
    iLink.edfdatafilename = strcat('ac_', int2str(participantNumber), '_vs');
    iLink = InitEyeLink(iLink, stimuliScrn);
end

%get filenames
trials= dir('im/search/*.png');
for t= 1:length(trials)
    trials(t).beta =  regexp(trials(t).name, '(?<=beta=)\d\.[\d]*', 'match');
    isTP = isempty(regexp(trials(t).name, 'TA', 'match'));
    trials(t).seed =  regexp(trials(t).name, '(?<=seed=)[\d]*', 'match');
    if isTP
        trials(t).isTP = 1;
        trials(t).x =  regexp(trials(t).name, '(?<=x=)[\d]*', 'match');
        trials(t).y =  regexp(trials(t).name, '(?<=y=)[\d]*', 'match');
        
    else
        trials(t).x = nan;
        trials(t).y = nan;
        trials(t).isTP = 0;
    end
end

% randomly permute
trials = trials(randperm(length(trials)));
pracTrials = trials(1:20);
trials = trials(randperm(length(trials)));

%% display intro screen
setupBlankAndFixCross;
Screen('DrawTexture', stimuliScrn, params.t_blank);
DrawFormattedText(stimuliScrn, 'First we will try 20 practise trials.', 'center', 'center')
Screen('Flip', stimuliScrn)
WaitSecs(0.5);
KbWait;
GetSecs; % run to load function into memory
for t= 1:length(pracTrials)
    rt = DoATrial(pracTrials(t), stimuliScrn, params, iLink);
end

Screen('DrawTexture', stimuliScrn, params.t_blank);
DrawFormattedText(stimuliScrn, 'Press any key when you are ready to start the experiment.', 'center', 'center')
Screen('Flip', stimuliScrn)
WaitSecs(0.5);
KbWait;


fout = fopen(strcat('results/ac_', int2str(participantNumber), '_vs.txt'), 'w');
fprintf(fout, 'participantNumber, trialNum, TargetPresent, beta, seed, x, y, RT\n');
for t= 1:length(trials)
    % check if it is time to give observer a break
    if t==81
        Screen('DrawTexture', stimuliScrn, params.t_blank);
        DrawFormattedText(stimuliScrn, 'You have completed 1 out of 3 blocks!', 'center', 'center')
        Screen('Flip', stimuliScrn)
        WaitSecs(0.5);
        KbWait;
    elseif t==160
        Screen('DrawTexture', stimuliScrn, params.t_blank);
        DrawFormattedText(stimuliScrn, 'Only one block left!', 'center', 'center')
        Screen('Flip', stimuliScrn)
        WaitSecs(0.5);
        KbWait;
    end
    rt = DoATrial(trials(t), stimuliScrn, params, iLink);
    if trials(t).isTP==1
        fprintf(fout, '%d, %d, %d, %s, %s, %s, %s, %.3f\n', participantNumber, t,...
            trials(t).isTP, trials(t).beta{1}, trials(t).seed{1}, trials(t).x{1}, trials(t).y{1}, rt);
    else
        fprintf(fout, '%d, %d, %d, %s, %s, NaN, NaN, %.3f\n', participantNumber, t, trials(t).isTP, trials(t).beta{1}, trials(t).seed{1}, rt);
    end
end
fclose(fout);
Screen('DrawTexture', stimuliScrn, params.t_blank);
DrawFormattedText(stimuliScrn, 'Thank You!', 'center', 'center')
Screen('Flip', stimuliScrn)


if iLink.doILink
    EYELINK('ReceiveFile',[iLink.edfdatafilename]);
    EYELINK('closefile');
    EYELINK('shutdown');
end

ShowCursor;
KbWait;

sca
end

function rt = DoATrial(trial, stimuliScrn, params, iLink)

% load image
im = imread(['im/search/' trial.name]);
t_im = Screen('MakeTexture', stimuliScrn, im);


% show fixation cross
Screen('DrawTexture', stimuliScrn, params.t_fixCross);
Screen('Flip', stimuliScrn)
% WaitSecs(params.fixationCrossDuration);
EyelinkDoDriftCorrection(iLink.el, params.midX, params.midY, 0, 1);
% Start data recording to EDF file, BEFORE params. */
status=Eyelink('startrecording');

% show stimulus
Screen('DrawTexture', stimuliScrn, t_im);
t0 = GetSecs;
Screen('Flip', stimuliScrn);

if iLink.doILink
    Eyelink('message', 'params_ON');	 % message for RT recording in analysis
    Eyelink('message','SYNCTIME');
end

if trial.isTP
    time2wait = params.maxStimuliparams;
else
    % decide on max stimulus params (shorter than usual as there is no
    % target!)
    if strcmp(trial.beta, '1.6')
        time2wait = params.catchTrialparams (1);
    elseif strcmp(trial.beta, '1.65')
        time2wait = params.catchTrialparams (2);
    elseif strcmp(trial.beta, '1.7')
        time2wait = params.catchTrialparams (3);
    end
end
resp = getObserverInput(time2wait);
rt = GetSecs-t0;
if iLink.doILink
    Eyelink('message','TRIAL_OVER');
end

%% if participant didn't find the target - blue screen of doom!
if (resp == 1) && (trial.isTP)
    Screen('DrawTexture', stimuliScrn, params.t_blank);
    Screen('Flip', stimuliScrn)
    WaitSecs(params.postTrialTime);
elseif  (resp == 0) && (trial.isTP)
    Screen('DrawTexture', stimuliScrn, params.t_blue);
    Screen('Flip', stimuliScrn)
    WaitSecs(params.penaltyTime)
    rt = -1;
elseif (resp == 1) && (~trial.isTP)
    Screen('DrawTexture', stimuliScrn, params.t_red);
    Screen('Flip', stimuliScrn)
    WaitSecs(params.penaltyTime)
elseif  (resp == 0) && (~trial.isTP)
    % show green screen
    Screen('DrawTexture', stimuliScrn, params.t_green);
    Screen('Flip', stimuliScrn)
    WaitSecs(params.postTrialTime);
    rt = -1;
end
Screen('Close', t_im);
end
