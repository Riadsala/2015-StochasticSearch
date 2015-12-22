function out = WaitForCentralFixation(iLink, display)

timeSinceStartedTrying = tic;
waitingForCentralFixation = true;

while waitingForCentralFixation
    timeWaiting = toc(timeSinceStartedTrying);
    [fx fy] = GetCurrentFixLoc(iLink);
    % work in dist^2, as that's a little faster
    d2 = (fx-display.midX)^2 + (fy-display.midY)^2;
    if d2 < iLink.centralTol
        % start timer - observer must fixate centre for x ms
        timeSinceStartedFixatingCentre = tic;
        while d2 < iLink.centralTol
            timeSpentFixatingCentre = toc(timeSinceStartedFixatingCentre);
            if timeSpentFixatingCentre > iLink.centralFixDur
                % we have a central fixation, so we can carry on!
                waitingForCentralFixation = false;
                out = 1;
                break; 
            end
            % get new distance to centre            
            [fx fy] = GetCurrentFixLoc(iLink);
            d2 = (fx-display.midX)^2 + (fy-display.midY)^2;
        end
        
    elseif timeWaiting > iLink.howLongToTryFor
        %% give up waiting and quit trial and recal
        waitingForCentralFixation = false;
        out = -1; % code to tell parent function to recalibrate!
    end
end
