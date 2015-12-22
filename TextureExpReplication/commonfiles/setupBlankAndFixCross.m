% blank screen
bkgrndGreyLevel = 128;
N = 512;
blank = bkgrndGreyLevel*ones(N);
params.t_blank = Screen('MakeTexture', stimuliScrn, blank);
blank(round(N/2-3):round(N/2+3),round(N/2-3):round(N/2+3)) = 64;
t_fixDot = Screen('MakeTexture', stimuliScrn, blank);
clear blank
% fixation cross
fixCross = makeFixCross(N, bkgrndGreyLevel);
params.t_fixCross = Screen('MakeTexture', stimuliScrn, fixCross);
clear fixCross
% penalty screen
red = 150*ones(N,N,3);
red(:,:,2:3) = 0;
params.t_red = Screen('MakeTexture', stimuliScrn, red);
clear red;
% reward screen
green = zeros(N,N,3);
green(:,:,2) = 150;
params.t_green = Screen('MakeTexture', stimuliScrn, green);
clear green;
