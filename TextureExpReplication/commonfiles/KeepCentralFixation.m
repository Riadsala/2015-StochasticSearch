function status = KeepCentralFixation(iLink, display, time2wait)

timer = tic;
if nargin < 3
    time2wait = iLink.stimuliDisplayTime;
end

while toc(timer) < time2wait
    [fx fy] = GetCurrentFixLoc(iLink);
    % work in dist^2, as that's a little faster
    d2 = (fx-display.midX)^2 + (fy-display.midY)^2;
    
    if d2 > iLink.centralTol
        status = false;
        return;
    end
    WaitSecs(0.01);
    
end

status = true;