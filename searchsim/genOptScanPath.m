function genOptScanPath
clear all
close all
N = 1024;
nFix = 50;
Beta = [1.60 1.65 1.70];
detModel = 'dir';
for b = 1:3
     scanpath{b} = runSim(N, Beta(b), detModel, nFix);
end
save optScanPath scanpath
end


function fix = runSim(N, beta, detModel, nFix)
 plotting = true;
init_fix = [N/2, N/2];

pTarg = ones(N,N);
X = PixelsToVisualAngle(repmat(1:N, [N 1]));
Y = X';

fix = [];
fix(1,:) = init_fix;

for t= 1:nFix
     tic
     disp(['working on fixation ' int2str(t)])
     % check how likely we were to miss the target, if it was there, at
     % every location
     Xt = X-PixelsToVisualAngle(fix(t,1));
     Yt = Y-PixelsToVisualAngle(fix(t,2));
     probOfMiss = reshape(1 - logisticReg_model(beta, Xt, Yt), [N,N]);
     pTarg = pTarg .* probOfMiss;
     clear Xt Yt probOfMiss
     
     % the best place to now fixated is the one that will minimise pTarg
     % we will consider a random sample of possible points
     F = randi(N, [1000 2]);
     for f = 1:length(F)
          Xt = X-PixelsToVisualAngle(F(f,1));
          Yt = Y-PixelsToVisualAngle(F(f,2));
          probOfMiss = reshape(1 - logisticReg_model(beta, Xt, Yt), [N,N]);
          pEffectOfFix(f) = sum(pTarg(:) .* probOfMiss(:));
          clear Xt Yt probOfMiss
     end
     
     [~,idx] = min(pEffectOfFix);
     fix(t+1,:) = F(idx,:); %#ok<AGROW>
     clear pEffectOfFix idx
     toc
end % trial
end



