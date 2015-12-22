

function fix = runSimRandom(N, beta, targLoc, alpha, maxFix, mkPlots)

% get saccade distribution to use
load ../clarke2009data/SaccDistribution.mat
Q = size(SaccByPos,1);

respArray = -0.5*ones(N, N);
respArray(targLoc(1), targLoc(2)) = 0.5;
init_fix = [N/2 N/2];

T=0;
targetfound = 0;
fix = [];
fix(1,:) = init_fix;


X = 0;

dE = 10;
prior = 1/N^2 * ones(N);

while targetfound==0

    T = T + 1;
    
    dprimeI(:,:,T) = calc_dprimeI(N, beta, fix(T,:));
    
    % generate location responses for this fixation
    W(:,:,T) = respArray + X + randn(N)*1./(dprimeI(:,:,T).^2);
    
    % now compute the posterior probability at each potential target location i after T fixations
   pTarg(:,:,T) = calc_probTarg(T, W, dE, dprimeI, prior);
      
    max(max(pTarg(:,:,T)));
    if max(max(pTarg(:,:,T))) > alpha
        %     if rand < probofdetection
        targetfound = 1;
        idx = find(pTarg(:,:,T) ==  max(max(pTarg(:,:,T))));
        [x, y] = ind2sub([N N], idx);
        fix(T+1,:) = [x y];
        
    elseif T==maxFix
        targetfound = 1;
        fix(T+1,:) = [NaN NaN];
    else
        % quantise current fixation
        F = ceil(fix(T,:)/4)
        F(F==0) = 1; % should be able to remove at a later date
        saccDist(:,:) = SaccByPos(F(1), F(2), :, :);
        C = cumsum(saccDist(:));
        r = sum(saccDist(:)) * rand;
        [~, idx] = min(abs(C-r));
        [y,x] = ind2sub([Q Q], idx);
        Fn = Q*max([x y], 1)- randi([0 31], [1 2]);
        fix(T+1,:) = round(Fn/8);
    end % trial
    if mkPlots
        im = imresize(pTarg(:,:,T),10, 'nearest');
        im = im./max(im(:));
        imshow(im,[])
        hold on

        plot(10*fix(:,1), 10*fix(:,2), '-rx', 'linewidth', 3, 'markersize', 15);        
         plot(10*targLoc(2), 10*targLoc(1), 'bo', 'markersize', 20);
        export_fig( ['fix_' int2str(T) '.png']);
       
        hold off
    end
    %     hold off
    %       im = imresize(W(:,:,T),10, 'nearest');
    %     im = im./max(im(:));
%     imshow(im,[])   
%     hold on
%      plot(10*fix(:,1), 10*fix(:,2), '-ro', 'linewidth', 3, 'markersize', 15); 
% 
%     export_fig( ['W_' int2str(T) '.png']);

end

if mkPlots
a = pTarg(targLoc(1), targLoc(2),:);
hold off
plot(a(:))
hold all
b = mean(mean(pTarg(:,:,:)));
plot(b(:))
end

end


function pT = calc_probTarg(T, W, dE, dI, prior)
g = 0;
for t = 1:T
    g = g + calc_g(dE, dI, T, t) .* W(:,:,t);
end
pT = prior .* exp(g);

pT = pT ./ sum(pT(:));
end


function g = calc_g(dE, dI, T, t)
fracT = dE.^2 * dI(:,:,t).^2;
fracB = dE.^2;
for x = 1:T
    fracB = fracB + dI(:,:,x).^2;
end
g = fracT ./ fracB;
end



function dprime = calc_dprimeI(N, b, fix)

v = repmat(1:N, [N, 1]);
x = (v - fix(1))./N;
y = (v' - fix(2))./N;
r = sqrt(x.^2 + y.^2);


rb = [-2.96066, -2.96066-2.49720,-2.96066-3.86623];
bb = [ 1.12441,  1.12441+1.30533,  1.12441+2.26095];
% Coefficients:
%              Estimate Std. Error t value Pr(>|t|)    
% (Intercept)   1.12441    0.09216  12.201  < 2e-16 ***
% r            -2.96066    0.37376  -7.921 2.08e-13 ***
% betamedium    1.30533    0.13033  10.015  < 2e-16 ***
% betasmooth    2.26095    0.13033  17.348  < 2e-16 ***
% r:betamedium -2.49720    0.52858  -4.724 4.54e-06 ***
% r:betasmooth -3.86623    0.52858  -7.314 7.47e-12 ***
dprime = bb(b) + rb(b).*r;
dprime(isinf(dprime(:))) = max(dprime(isfinite(dprime)));
end
