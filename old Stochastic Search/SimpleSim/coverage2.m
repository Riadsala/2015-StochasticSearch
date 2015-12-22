function mfi = CoverageMetricFunction(coords, N)

numfix = min(size(coords,1),50);

coords = coords/4;
N=N./4;

dists = zeros(N,N,numfix);
for t=1:numfix
    X=4*repmat((1:N)-coords(t,1), [N,1]);
    Y=4*repmat((1:N)'-coords(t,2), [1,N]);
    dists(:,:,t) = sqrt(X.^2+Y.^2);    
end
for t=1:numfix
    mindists(:,:,t) = dists(:,:,1);
    for n=1:t
        mindists(:,:,t) = min(mindists(:,:,t), dists(:,:,n));
        mfi(t) = mean(mean(mindists(:,:,t)));
    end
end

mfi(1:(numfix-1)) = mfi(1:(numfix-1)) - mfi(2:(numfix));
mfi(numfix)=[];
% plot(mfi)
% xlabel('fixation number');
% ylabel('mean distance from pixel to nearest fixation');
% title('Independant Random Searcher');
