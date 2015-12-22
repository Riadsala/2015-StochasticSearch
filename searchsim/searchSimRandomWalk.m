function runExperiments
clear all
close all
N = 1024;
reps = 5;
R = [100 225 350];
Phi = 0:45:315;
Beta = [1.60 1.65 1.70];
detModel = 'iso';

for i = 1:reps
    for j = 1:length(R)
        for k = 1:length(Phi)
            for b = 1:3

                r = R(j);
                phi = Phi(k);
                beta = Beta(b);
                res.ran.F{i,j,k,b} = runSimRandom(N, beta, r, phi, detModel);
                res.ran.nFix(i,j,k,b) = size(res.ran.F{i,j,k,b},1);
                res.opt.F{i,j,k,b} = runSimOptimal(N, beta, r, phi, detModel);
                res.opt.nFix(i,j,k,b) = size(res.opt.F{i,j,k,b},1);
            end
        end
    end
end 

save searchSimResults res
resultsToText;
end


function fix = runSimRandom(N, beta, r, phi, detModel)

% get saccade distribution to use
load ../clarke2009data/SaccDistribution.mat
Q = size(SaccByPos,1);
% place target
[x y] = pol2cart(phi, r);
x = ceil(x); y = ceil(y);
targetloc = [x y]+N/2;

init_fix = [N/2, N/2];
% add a slight random offset, so there's even chance of us fixating any of
% the four centre points
init_fix = init_fix + randi(2,[1 2])-1;

t=0;
targetfound = 0;
fix = [];
fix(1,:) = init_fix;

while targetfound == 0
 t = t+1;
    dist2targetX = PixelsToVisualAngle((fix(t,1)-targetloc(1)));
    dist2targetY = PixelsToVisualAngle((fix(t,2)-targetloc(2)));
    
    probofdetection = logisticReg_model(beta, dist2targetX, dist2targetY);
    clear dist2targetX dist2targetY
    

    if rand < probofdetection
        targetfound = 1;
        fix(t+1,:) = targetloc; %#ok<AGROW>
    elseif t==50
        targetfound = 1;
        fix(t+1,:) = [NaN NaN];
    else
        % quantise current fixation
        F = ceil(fix(t,:)/Q);
        saccDist(:,:) = SaccByPos(F(1), F(2), :, :);
        C = cumsum(saccDist(:));
        r = sum(saccDist(:)) * rand;
        [~, idx] = min(abs(C-r));
        [x y] = ind2sub([Q Q], idx);
        Fn = Q*max([x y], 1)- randi([0 31], [1 2]);;
        fix(t+1,:) = Fn;
    end % trial
%     plot(fix(:,1), fix(:,2), '-bx'); axis([1 1024 1 1024])
%     hold on
%     plot(targetloc(1), targetloc(2), 'rx', 'markersize', 10);
end

end



function fix = runSimOptimal(N, beta, r, phi, detModel)

% get saccade distribution to use
load optScanPath

% place target
[x y] = pol2cart(phi, r);
x = ceil(x); y = ceil(y);
targetloc = [x y]+N/2;

init_fix = [N/2, N/2];
% add a slight random offset, so there's even chance of us fixating any of
% the four centre points
init_fix = init_fix + randi(2,[1 2])-1;

t=0;
targetfound = 0;
fix = [];
fix(1,:) = init_fix;
if beta==1.6
    b = 1;
elseif beta == 1.65
    b = 2;
elseif beta == 1.7
    b = 3;
end
    
while targetfound == 0
 t = t+1;
    dist2targetX = PixelsToVisualAngle((fix(t,1)-targetloc(1)));
    dist2targetY = PixelsToVisualAngle((fix(t,2)-targetloc(2)));
    
    probofdetection = logisticReg_model(beta, dist2targetX, dist2targetY);
    clear dist2targetX dist2targetY
    

    if rand < probofdetection
        targetfound = 1;
        fix(t+1,:) = targetloc; %#ok<AGROW>
    elseif t==50
        targetfound = 1;
        fix(t+1,:) = [NaN NaN];
    else
        % just take next saccade from precomputed optimal
        fix(t+1,:) = scanpath{b}(t+1,:);
    end % trial

end

end


