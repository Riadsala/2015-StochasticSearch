
function [resp, responseKeyHit] = getObserverInputFS

% check if a key has been pressed
[keyIsDown, ~, keyCode] = KbCheck;
% wait for keypress
while ~keyIsDown
     [keyIsDown, ~, keyCode] = KbCheck;
end
% check what key was pressed
if find(keyCode) == KbName('j');
     responseKeyHit = 1;
     resp = 1;
elseif find(keyCode) == KbName('f');
     responseKeyHit = 1;
     resp = 0;
elseif  find(keyCode) == KbName('q');
     responseKeyHit = 1;
     resp = 3;
else
     resp = 3;
     responseKeyHit = 0;
end
end
