function [t results] = LookAtResp(resp, iLink, t, results, x, y, clp, targPresent)
if resp == 3
    sca;
    ShowCursor;
    return;
elseif resp == -1
    % failed to get stable central fixation - recalibrate
    disp('getting new calib');
    EyelinkDoTrackerSetup(iLink.el);
elseif resp == -2
    % user broke central fixation - discard trial
    t = t-1;
else
    if nargin < 4
        results = [];
    elseif nargin == 8
        results(t, :) = [targPresent x y a resp clp];
    elseif nargin == 4
        results(t, :) = resp;
    end
end
end