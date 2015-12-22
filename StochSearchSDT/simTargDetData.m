
clear all
close all


N = 1024;


alpha = 0.035; % stopping criteria

for beta=1:3
    j=0;
for i = 0:64:N
    j = j+1;
    targLoc = [i N/2];  
    for r = 1:100
        fix = runSimRandom2(N, beta, targLoc, alpha, 1, false);
        if sum(fix(2,:) == targLoc)==2
            res(beta,j,r) = 1;
            fp(beta,j,r) = 0;
        elseif isfinite(fix(2,1))
            fp(beta,j,r) = 1;
            res(beta,j,r) = 0;
        else
            res(beta,j,r) = 0;
            fp(beta,j,r) = 0;
        end
    end
end
end

save searchSimResults res
% resultsToText;
plot(([0:64:N]-512)*0.0139, mean(res(:,:,:),3)')
xlabel('eccentricity (visual degrees)');
xlim([-6 6])
ylabel('accuracy')
hold all
