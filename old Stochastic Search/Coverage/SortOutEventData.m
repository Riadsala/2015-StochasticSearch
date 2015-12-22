function [Results timeoffsets] = SortOutEventData(Person)
timeoffset = 0;
timeoffsets = [0 0 0];
% Results = zeros(300,13);
for erun = 1:3
    timeoffsets(erun)= timeoffset;
    filename = [Person int2str(erun) 'EVD.txt'];
    EventFile = fopen(filename);
    data = fscanf(EventFile, '%c',inf);
    fclose(EventFile);
    %%% delete header
    cursor = regexp(data, 'Time', 'once');
    data(1:(cursor-1)) = [];
    N = size(data,2);
    % find fixcross.jpg.
    FixcrossLocations = regexp(data, 'fixcross.jpg');
    n = size(FixcrossLocations, 2);
    pngLocations = regexp(data, 'png');
    % for each fixcross
    % check time for following imageshow event (as it could be a rouge
    % keypress!)
    % then get start and stop time for trial.
    for t = 1:n
        cursor = FixcrossLocations(t);
        cursor = cursor + 13;
        followingevent = regexp(data(cursor:N), 'Keyboard|ShowSlide', 'once', 'match');
        % Catch any keypresses made while the Fixation Cross is
        % being displayed
        while strcmp(followingevent, 'Keyboard')
            cursor=regexp(data(cursor:N), 'Keyboard', 'once')+cursor;
            cursor=regexp(data(cursor+8:N), '[a-zA-Z]', 'once')+cursor;
            followingevent = regexp(data(cursor:N), 'Keyboard|ShowSlide', 'once', 'match');
        end
        if strcmp(followingevent, 'Keyboard')
            disp('bug to fix');
            break
        end
        StartTime = str2double(regexp(data(cursor:N), '[0-9][0-9][0-9][0-9]+', 'once', 'match'));
        StartTime = StartTime+timeoffset;
        Filename = regexp(data(cursor:N), 'RMS[._=\w]*png',  'once', 'match');
        cursor = pngLocations(t)+3;
        EndTime = str2double(regexp(data(cursor:N), '[0-9]*', 'once', 'match'));
        EndTime = EndTime+timeoffset;
        ReactionTime = (EndTime-StartTime)/1000;
        followingevent = regexp(data(cursor:N), 'Keyboard|ShowSlide', 'once', 'match');
        if strcmp(followingevent, 'ShowSlide')
            disp('bug to fix');
            break
        end
        Keypress =   regexp(data(cursor:N), '(?<=\s)[^0-9](?=\s)', 'once', 'match');
        Reply = size(regexp(Keypress, '[kyuiophjl;nm,.\[YUIOPHJKLNM]'),1)+1;
        Target = size(regexp(Filename, 'TP'),1)+1;
        Correct = (Reply==Target);
        RMS = str2double(regexp(Filename, '(?<=RMS=)(1(\.2)*)|0\.8(?=_)', 'match'));
        beta = str2double(regexp(Filename, '(?<=beta=)1\.[567]*(?=_seed)', 'match'));
        seed = str2double(regexp(Filename, '(?<=seed=)[0-9]*(?=_T)', 'match'));
        % Find target location in a suprisingly dumbway. Shrugs shoulders
     
        Results(t+(erun-1)*100).basicinfo(1:8) = [StartTime, EndTime, Target, RMS, beta, seed, Correct, ReactionTime];
    end
    timeoffset = timeoffset + EndTime;
end
end