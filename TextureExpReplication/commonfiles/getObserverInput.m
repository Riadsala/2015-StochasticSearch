function resp = getObserverInput(time2wait)
[keyIsDown, ~, keyCode] = KbCheck;
t0 = GetSecs;
% wait for keypress
responseKeyHit = 0;
while ~responseKeyHit
    while ~keyIsDown
        [keyIsDown, ~, keyCode] = KbCheck;
        WaitSecs(0.01);
        if GetSecs-t0 > time2wait
            resp = 0;
            responseKeyHit = 1;
            break;
        end
    end
    % check what key was pressed
    if find(keyCode) == KbName('space');
        responseKeyHit = 1;
        resp = 1;
    elseif find(keyCode) == KbName('f');
        responseKeyHit = 1;
        resp = 0;
    elseif  find(keyCode) == KbName('q');
        responseKeyHit = 1;
        resp = 3;
    end
end