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
q = 4;
pTarg = ones(N/q,N/q);
X = PixelsToVisualAngle(repmat(1:N/q, [N/q 1]));
Y = X';

fix = [];
fix(1,:) = init_fix;

for t= 1:nFix
     tic
     disp(['working on fixation ' int2str(t)])
     % check how likely we were to miss the target, if it was there, at
     % every location
     Xt = q*X-PixelsToVisualAngle(fix(t,2));
     Yt = q*Y-PixelsToVisualAngle(fix(t,1));
     probOfMiss = 1 - logisticReg_model(beta, Xt, Yt);
     pTarg = pTarg .* probOfMiss;
     clear Xt Yt probOfMiss
     
     % the best place to now fixated is the one that will minimise pTarg
     % we will consider a random sample of possible points
     
     for  x = q:q:1024
          for y = q:q:1024
               Xt = q*X-PixelsToVisualAngle(y);
               Yt = q*Y-PixelsToVisualAngle(x);
               probOfMiss =1-logisticReg_model(beta, Xt, Yt);
               updatedPotentialPTarg = pTarg.* probOfMiss;
               pEffectOfFix(x/q, y/q) = mean(updatedPotentialPTarg(:));
               clear Xt Yt probOfMiss
          end
     end
     [~,idx] = min(pEffectOfFix(:));
     [a b]= ind2sub([1024/q, 1024/q], idx);
     fix(t+1,:) = [q*a q*b]; %#ok<AGROW>
     clear a b idx
          
     if plotting
          subplot(1,2,1)
          contour(pTarg, 0:0.01:1)
          hold all          
          plot(fix(:,2)/q, fix(:,1)/q, 'ok-', 'linewidth', 2)          
          axis square
          title('current prob of target');
          subplot(1,2,2)
          pEffectOfFix(pEffectOfFix>quantile(pEffectOfFix(:), 0.75)) = NaN;
          contour(pEffectOfFix,100)
          
          hold all
          plot([0, 1024/q], [512/q, 512/q], 'k-')
          plot([512/q, 512/q],[0, 1024/q], 'k-')
          plot(fix(:,2)/q, fix(:,1)/q, 'ok-', 'linewidth', 2)
          axis square
          export_fig(['beta_' num2str(beta) '_fixation_' int2str(t) '.png'])
          close all
     end
     
     clear pEffectOfFix idx
     toc
end % trial
end

function p = logisticReg_model(beta, x, y)
r = sqrt(x.^2 + y.^2);
x2 =  x.^2;

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


% function p = linear_model(beta, r)
% p = -597 + 409.145*beta - 11.441*r;
% end
%
% function p = gauss_model(beta, x,y, type)
% if strcmp(type, 'iso')
%      r = sqrt(x.^2 + y.^2);
%      switch beta
%           case 1.6
%                a = 0.95;
%                z = 1.36316;
%                k = 0.70800;
%           case 1.65
%                a = 0.80191;
%                z = 4.02467;
%                k = 1.65307;
%           case 1.7
%                a = 0.95727;
%                z = 4.02467;
%                k = 1.94994;
%      end
%      p = a * exp(-(r./z).^k);
% else
%
%      switch beta
%           case 1.6
%                a = 0.485;
%                zx = 2.70;
%                zy = 3.455;
%           case 1.65
%                a = 0.76;
%                zx = 3.378;
%                zy = 5.23;
%           case 1.7
%                a = 0.9567;
%                zx = 4.397;
%                zy = 5.63;
%      end
%      p = a * exp(-((x./zx).^2+(y./zy).^2));
% end
% end
