function Coverage = CoverageMetricFunction(coords, N)


N=256;
numfix = 3;%size(coords,1);

SaccadeAmps = cumsum(csvread('SaccadeAmps.csv'));
SaccadeAmps(15:34) = [];
fix = [N/2, N/2];
coords=fix;
for n=1:numfix-1
    p = rand*max(SaccadeAmps);
    saccampmin = VisualAngleToPixels(min(find(SaccadeAmps>p))/2-0.5);
    saccampmax = VisualAngleToPixels(min(find(SaccadeAmps>p))/2);
    X = repmat((1:N)-fix(1), [N,1]);
    Y = repmat((1:N)'-fix(2), [1,N]);
    dists = sqrt(X.^2+Y.^2);
    dists = (dists<saccampmax).*(dists>saccampmin);
    if sum(dists(:))==0
        disp('oops')
    end
    p = sum(dists(:))*rand;
    dists = cumsum(dists(:));
    [fix(1) fix(2)] = ind2sub(N,find(dists<p, 1, 'last' ));
    coords = [coords; fix];
end


% figure; plot(coords(:,2), coords(:,1), ':rx');
% axis([1, N, 1,N]);
% hold on
area = zeros(N,N,numfix);
for t=1:numfix
    t
    for i = 1:N
        [i t]
        for j = 1:N
  
            for n=1:t
                disttofix(n) = norm([i,j] - coords(n,:));
            end
            area(i,j,t)=PixelsToVisualAngle(min(disttofix));
        end
    end
end

Coverage = mean(area(:));
contour(area, [1,2,3,4,5])
figure
hist(area(:))
