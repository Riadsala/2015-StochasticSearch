

function fix = runSimRandom(N, beta, targLoc, alpha, maxFix, mkPlots)


% get saccade distribution to use
load ../clarke2009data/SaccDistribution.mat
Q = size(SaccByPos,1);

x = meshgrid(0:64:N);
y = x';
respArray = [x(:), y(:)];
clear x y
% assign target to closet fit
d = (respArray(:,1)-targLoc(1)).^2 + (respArray(:,2)-targLoc(2)).^2;
[~, idx] = min(d);
clear d;
respArray(:,3) = -0.5;
respArray(idx,3) = 0.5;
clear idx


init_fix = [N/2 N/2];

T=0;
targetfound = 0;
fix = [];
fix(1,:) = init_fix;


X = 0;

dE = 10;
prior = 1/length(respArray) * ones(1, length(respArray));

while targetfound==0
    
    T = T + 1;
    
    dprimeI(:,T) = calc_dprimeI(respArray, beta, fix(T,:));
    
    % generate location responses for this fixation
    W(:,T) = respArray(:,3) + X + randn(length(dprimeI),1)*1./(dprimeI(:,T).^2);
    
    % now compute the posterior probability at each potential target location i after T fixations
    pTarg(:,T) = calc_probTarg(T, W, dE, dprimeI, prior);
    
    max(max(pTarg(:,T)));
    if max(pTarg(:,T)) > alpha
        %     if rand < probofdetection
        targetfound = 1;
        idx = find(pTarg(:,T) ==  max(pTarg(:,T)));
         [respArray(idx,1), respArray(idx,2)] ;
        fix(T+1,:) = [respArray(idx,1), respArray(idx,2)] ;
        
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
    g = g + calc_g(dE, dI, T, t) .* W(:,t);
end
pT = prior .* exp(g');

pT = pT ./ sum(pT(:));
end


function g = calc_g(dE, dI, T, t)
fracT = dE.^2 * dI(:,t).^2;
fracB = dE.^2;
for x = 1:T
    fracB = fracB + dI(:,x).^2;
end
g = fracT ./ fracB;
end



function dprime = calc_dprimeI(respArray, b, fix)

d =  sqrt((respArray(:,1)-fix(1)).^2 + (respArray(:,2)-fix(2)).^2);

bb = [1.052e+00 1.052e+00+8.314e-01 1.052e+00+1.506e+00];
xb = [-7.073e-06 -7.073e-06-1.002e-06 -7.073e-06-2.559e-06];
yb = -1.085e-05;
% (Intercept)    1.052e+00  9.163e-02  11.479  < 2e-16 ***
% x2            -7.073e-06  1.526e-06  -4.636 9.59e-06 ***
% betamedium     8.314e-01  1.175e-01   7.074 1.34e-10 ***
% betasmooth     1.506e+00  1.175e-01  12.815  < 2e-16 ***
% y2            -1.085e-05  9.138e-07 -11.876  < 2e-16 ***
% x2:betamedium -1.002e-06  2.116e-06  -0.474    0.637    
% x2:betasmooth -2.559e-06  2.116e-06  -1.209    0.229    
dprime = bb(b) + xb(b)*(respArray(:,1)-fix(1)).^2 + yb*(respArray(:,2)-fix(1)).^2;
dprime(dprime<0.1) = 0.1;
end
