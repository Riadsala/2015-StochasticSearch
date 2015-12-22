function EstimateSaccDist
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% .m file to read, process, and display results from Eye tracking
% experiment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saccDat = csvread('saccDat.txt', 1);

Q = 32;
SaccByPos = zeros(Q, Q, Q, Q);


sacc = ceil(saccDat(:,5:8)./Q);

% don't use last saccade, as that's probably towards the target
for s = 1:length(saccDat)
    SaccByPos(sacc(s,1),sacc(s,2),sacc(s,3),sacc(s,4)) = ...
        SaccByPos(sacc(s,1),sacc(s,2),sacc(s,3),sacc(s,4)) + 1;
    
end


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
