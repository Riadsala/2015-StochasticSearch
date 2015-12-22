function numrefixes = checkForReFixations(X)


% distance for refixation
fix_radius = 60;
% how far in the past to check for re-fixations


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numrefixes = 0;
numfix = size(X,1);
stepstocheck = 3;
% for all fixations after the first...
for t=2:(numfix-2)
    %  go back through the last n fixations to check for refixations
    for f = max((t-stepstocheck),1):(t-1)
        interfixdistance = norm(X(t,:)-X(f,:));
        numrefixes = numrefixes + (interfixdistance<fix_radius);
    end
end
