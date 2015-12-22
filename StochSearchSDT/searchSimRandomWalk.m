
clear all
close all


N = 128;


alpha = 0.99; % stopping criteria
targLoc = [16 64];
for i = 1:100
fix = runSimRandom(N, 1, targLoc, alpha, 50, true);
res(i) = length(fix);
end

boxplot(res)
% resultsToText;

