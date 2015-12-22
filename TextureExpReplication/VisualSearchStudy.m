function VisualSearchStudy
addpath('commonfiles/');
participantNumber = input('enter participant number:');

iLink.doILink = 0;

%% set things up
bkgrndGreyLevel = 127;
N = 512;
[stimuliScrn] = Screen('OpenWindow',0, bkgrndGreyLevel);%, [001 01 1600 900]
[display.width, display.height]=Screen('WindowSize', stimuliScrn);%screen returns the size of the window
display.midX = round(display.width/2);
display.midY = round(display.height/2);
% HideCursor;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize all eyelink and calibrate
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iLink.doILink == 1
     iLink = InitEyeLink(iLink, stimuliScrn);
     iLink.howLongToTryFor = 3;
     iLink.centralFixDur = fixCrossDisplayTime;
     iLink.stimuliDisplayTime = stimuliDisplayTime;
end

%get filenames
trials= dir('im/search/*.png');
for t= 1:604
trials(t).beta =  regexp(trials(t).name, '(?<=beta=)\d\.[\d]*', 'match');
isTP = isempty(regexp(trials(t).name, 'TA', 'match'));
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


foutWW = fopen('results.txt', 'w');
fprintf(foutWW, 'participantNumber trialNum, TargetPresent, x, y RT\n');
for t= 1:5
     rt = DoATrial(trials(t), stimuliScrn);
     fprintf(foutWW, '%d, %d, %d, %d, %d, %d\n', participantNumber, t, trials(t).isTP, trials(t).x, trials(t).y, rt);
end
fclose(foutWW);


sca
end

function click = DoATrial(trial, stimuliScrn)
fixCrossDuration = 1;

% blank screen
bkgrndGreyLevel = 127;
N = 512;
blank = bkgrndGreyLevel*ones(N);
t_blank = Screen('MakeTexture', stimuliScrn, blank);
blank(round(N/2-3):round(N/2+3),round(N/2-3):round(N/2+3)) = 64;
t_fixDot = Screen('MakeTexture', stimuliScrn, blank);
clear blank
% fixation cross
fixCross = makeFixCross(N, bkgrndGreyLevel);
t_fixCross = Screen('MakeTexture', stimuliScrn, fixCross);
clear fixCross

% show fixation cross
Screen('DrawTexture', stimuliScrn, t_fixCross);
Screen('Flip', stimuliScrn)
% load image

im = imread(['im/search/' trial.name]);
t_im = Screen('MakeTexture', stimuliScrn, im);
Screen('DrawTexture', stimuliScrn, t_im);

WaitSecs(fixCrossDuration);
tic
Screen('Flip', stimuliScrn)
resp = getObserverInput;
click = toc;
Screen('DrawTexture', stimuliScrn, t_blank);
Screen('Flip', stimuliScrn)
WaitSecs(0.5);



end


