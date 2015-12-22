function runExperiments
clear all
close all
N = 1024;
reps = 1;
R = [100 225 350];
Phi = 0:45:315;
Beta = [1.60];
detModel = 'iso';

for i = 1:reps
    for j = 1:length(R)
        for k = 1:length(Phi)
            for b = 1
                [i j k b]
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

% save
fout = fopen('exScanpathsForPlot.txt', 'w');
fprintf(fout, 'sim, t, x, y\n')
for t = 1:20
     fprintf(fout, 'optimal, %d, %d, %d\n', t, res.opt.F{1}(t,1),res.opt.F{1}(t,2))
end
for t = 1:20
     fprintf(fout, 'stochastic, %d, %d, %d\n', t, res.ran.F{1}(t,1),res.ran.F{1}(t,2))
end
fclose(fout)
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

for t = 1:20
    dist2targetX = PixelsToVisualAngle((fix(t,1)-targetloc(1)));
    dist2targetY = PixelsToVisualAngle((fix(t,2)-targetloc(2)));
    
    probofdetection = logisticReg_model(beta, dist2targetX, dist2targetY);
    clear dist2targetX dist2targetY
    


        % quantise current fixation
        F = ceil(fix(t,:)/Q);
        saccDist(:,:) = SaccByPos(F(1), F(2), :, :);
        C = cumsum(saccDist(:));
        r = sum(saccDist(:)) * rand;
        [~, idx] = min(abs(C-r));
        [x y] = ind2sub([Q Q], idx);
        Fn = Q*max([x y], 1)- randi([0 31], [1 2]);;
        fix(t+1,:) = Fn;
  
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
    
for t= 1:20
    dist2targetX = PixelsToVisualAngle((fix(t,1)-targetloc(1)));
    dist2targetY = PixelsToVisualAngle((fix(t,2)-targetloc(2)));
    
    probofdetection = logisticReg_model(beta, dist2targetX, dist2targetY);
    clear dist2targetX dist2targetY
    

 
        % just take next saccade from precomputed optimal
        fix(t+1,:) = scanpath{b}(t+1,:);
  

end

end



function p = logisticReg_model(beta, x, y)
r = sqrt(x.^2 + y.^2);
x2 = x.^2;

switch beta
     case 1.6
          b1 = 0.158655;
          br = -0.500250;
          bx = 0.003545;
     case 1.65
          b1 =  1.617402+0.158655;
          br = -0.500250 -0.231463;
          bx = 0.003545+0.028478;
     case 1.7
          b1 =  2.760698+0.158655;
          br = -0.500250 -0.247090;
          bx = 0.003545+ 0.015953;
end
 lm = b1 + br.*r + bx.*x2;

p = 1./(1+exp(-lm));

end

