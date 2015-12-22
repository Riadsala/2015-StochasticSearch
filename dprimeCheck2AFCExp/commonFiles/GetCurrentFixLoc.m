function [x y output] = GetCurrentFixLoc(iLink)
%%%%%%%%%%%%%%%%%%%%%%%%
%  grab current eye position
if iLink.doILink
    if EYELINK('isconnected')==iLink.el.notconnected   % Check link often so we don't lock up if tracker lost
        output.message = 'ForceQuitProgram';
        x = NaN; y = NaN;
        return;
    end;
    
    %ensure we maintain fixation.  Just check location
    while Eyelink( 'NewFloatSampleAvailable') < 1        
        WaitSecs(0.01);
    end
    evt = Eyelink( 'NewestFloatSample');
    x = (max(evt.gx));
    y = max(evt.gy);
 
    
end