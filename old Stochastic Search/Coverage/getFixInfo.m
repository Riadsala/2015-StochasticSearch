 function  Fixations = getFixInfo(StartTime, EndTime, Fixs)
        FixDurationParam = 100;
        FinalFixMinDuration = 200;
        
xoffset = (1200-1024)/2;
yoffset = (1600-1024)/2;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Note: I am including the final fixation on the fixation cross as
        % a trial data point... Otherwise I have no data for my first saccade!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        index= find((Fixs(:,2)<=EndTime).*(Fixs(:,2)>=StartTime));
        if size(index,1)>0
            % sort out initial fixation, which probably started
            % before the onset of the stimulus. We will consider the time
            % from stimulus onset to end of fixation as a fixation. Unless
            % it lasts shorter than 200ms.
            if index(1)>1
                index = [index(1)-1; index];
            end
            CroppedInitFixLength = Fixs(index(1),2)+Fixs(index(1),3)-StartTime;
            if (CroppedInitFixLength<FixDurationParam)*(size(index,1)>1)
                index(1) = [];
            end
            % now sort out final fixation. Only count if participant
            % fixates for at least 200ms before hitting button.
            TimeBetweenFinalFixAndButtonPress = EndTime- Fixs(index(size(index,1)),2);
            if  (TimeBetweenFinalFixAndButtonPress < FinalFixMinDuration)*(size(index,1)>1)
                index(size(index,1)) = [];
            end
            Fixations.y = 1024-(1200-Fixs(index, 5)-xoffset);
            Fixations.x = Fixs(index, 4)-yoffset;
%             Fixations.coverage = coverage([Fixations.x, Fixations.y],1024);           
            Fixations.number= size(index,1);
    
        else
            % no fixations recorded for trial!!!
            Fixations.number = 0;
            Fixations.coverage = 0;
            Fixations.x = [];
            Fixations.y = [];
        end
    end % Fixations