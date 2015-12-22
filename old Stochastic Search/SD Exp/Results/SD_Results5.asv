function SD_Results4T
close all
clear all
% %define offset to turn Clearview output coords to image coords.
xoffset = (1200-1024)/2;
yoffset = (1600-1024)/2;
BlocksInRun = 33;
TrialsInBlock = 4;
TrialCounter = 0;
% Results = zeros(2160, 5);
set(0,'DefaultAxesColorOrder',[0,0,0])
set(0,'DefaultAxesLineStyleOrder',{'-o','->','-s', '--o','-->','--s'})
Person = {'LM', 'MY'};
preceedTrialStartTime = 50;
postTrialEndTime = -100;
maxAllowedDist=60;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sort out EventData - filenames and keypresses
%%%%%%help %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drullert = 0;
figure('Position',[1 1 400 800])
for Pctr=1:2
    TrialCounter = 0;
    Results = [];
    for part = 1:20
        FixFile = fopen( [Person{Pctr} int2str(part) 'EFD.txt']);
        FixData = fscanf(FixFile, '%c',inf);
        fclose(FixFile);
        %%% delete header
        cursor = regexp(FixData, 'Timestamp', 'once');
        FixData(1:(cursor-1)) = [];
        cursor = regexp(FixData, '\d', 'once');
        FixData(1:(cursor-1)) = [];
        % convert to matrix - first need to replace validity strings with
        % numbers
        FixData=regexprep(FixData, 'Both', '1');
        FixData=regexprep(FixData, 'None', '0');
        FixData=regexprep(FixData, 'LeftOnly', '0');
        FixData=regexprep(FixData, 'RightOnly', '0');
        Fixs = str2num(FixData);
        clear FixData
        EventFile = fopen([Person{Pctr} int2str(part) 'EVD.txt']);
        EventData = fscanf(EventFile, '%c',inf);
        fclose(EventFile);
        %%% delete header
        cursor = regexp(EventData, 'Description', 'once');
        EventData(1:(cursor+14))=[];
        N = size(EventData,2);
        % find fixcross.jpg.
        FixcrossLocations = regexp(EventData, 'fixcross.jpg');
        %     n = size(FixcrossLocations, 2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % go through all the trials in the run and get break times
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ImageLocations = regexp(EventData, '_beta');
        T = size(ImageLocations,2);
        for b = 1:BlocksInRun
            for t = 1:TrialsInBlock
                TrialCounter = TrialCounter+1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % find relvant fixation cross and get following image display time
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                cursor = FixcrossLocations((b-1).*TrialsInBlock+t);
                TrialStartTime = str2double(regexp(EventData(cursor:N), '\s\d+\s', 'match', 'once'))-preceedTrialStartTime;
                cursor = ImageLocations((b-1).*TrialsInBlock+t);
                TrialEndTime = str2double(regexp(EventData(cursor:N), '\s\d+\s', 'match', 'once'))+postTrialEndTime;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Check fixation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Fidx = find((Fixs(:,1)<TrialEndTime).*(Fixs(:,1)>TrialStartTime));
                Fixations.y = 1024-(1200-Fixs(Fidx, 4)-xoffset);
                Fixations.x = Fixs(Fidx, 3)-yoffset;
                DistancesFromCentre = getDistances(Fixations);
                Results(TrialCounter, 5) = mean(DistancesFromCentre);
                Results(TrialCounter, 6) = std(Fixations.x);
                Results(TrialCounter, 7) = std(Fixations.y);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get surface parameters from filename
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                cursor = ImageLocations((b-1).*TrialsInBlock+t);
                Filename = regexp(EventData(cursor:N), 'beta[._=\w]*png',  'once', 'match');
                beta = str2double(regexp(Filename, '(?<=beta=)1\.[567]*', 'match'));
                r = str2double(regexp(Filename, '(?<=r=)(50|100|150|200|250|300|350|400|450)', 'match'));
                if r>0
                    phi = str2double(regexp(Filename, '(?<=phi=)(0|45|90|135|180|225|270|315)', 'match'));
                else
                    r = inf;
                    phi = inf;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % new filter - include trials with mean,std fix < threshold
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                dRule = (Results(TrialCounter, 5)<60)*(Results(TrialCounter, 6)<40)*(Results(TrialCounter, 7)<40);
                drullert = drullert +  dRule;
                if dRule&(size(Fidx,1)>0)
                    % get Right|Left Keyboard response
                    Response =  regexp(EventData(cursor:N), '(Right)|(Left)',  'once', 'match');
                    % check that response was found before the next trial was shown
                    % if not, set response to 2 (do not count)
                    if (b-1)*TrialsInBlock+t~=T
                        RespLoc =  regexp(EventData(cursor:N), '(Right)|(Left)',  'once');
                        if (RespLoc+cursor)>ImageLocations((b-1).*TrialsInBlock+t+1)
                            %                         disp('problem with missing response');
                            %                         disp(part)
                            %                         disp(Filename)
                            Results(TrialCounter,4) = 2;
                        else
                            if strcmp(Response, 'Right')
                                Results(TrialCounter,4) = 0;
                            elseif strcmp(Response, 'Left')
                                Results(TrialCounter,4) = 1;
                            end
                        end
                    end
                elseif size(Fidx,1)==0
                    Results(TrialCounter,4) = 3;
                else
                    Results(TrialCounter,4) = 2;
                end
                Results(TrialCounter,1:3) = [beta, r, phi];
            end
        end
    end
    Pctr    
    subplot(2,1,Pctr)
    betaSet = [1.6, 1.65, 1.7];
    rSet = 50:50:450;
    rSet = rSet;
    phi = [0,45,90,135,180,225,270,315];
    %     Accuracy = zeros(3,8,8);
    AccuracyOverEcc = zeros(3,9);
    TargetAbsentAccuracy = zeros(3,1);
    TA2Plot = zeros(3,9);
    for b = 1:3
        for r=1:9
            for p=1:8
                beta = betaSet(b);
                index = find((Results(:,1)==beta).*(Results(:,2)==rSet(r)).*(Results(:,3)==phi(p)).*(Results(:,4)~=2));
                Accuracy(b,r,p, Pctr) = sum(Results(index, 4)==1)./size(Results(index, 4)==1, 1);
            end
            tmp = Accuracy(b,r,:,Pctr);
            tmp=tmp(:);
            AccuracyOverEcc(b,r) = mean(tmp(isfinite(tmp)));
        end
        index = find((Results(:,1)==beta).*isinf(Results(:,2)).*(Results(:,4)~=2));
        TargetAbsentAccuracy(b) = 1-mean(Results(index, 4)); %#ok<FNDSB>
        TA2Plot(b,:) = TargetAbsentAccuracy(b)*ones(1,9);
      
        PersonAcc(:,:,Pctr) =  AccuracyOverEcc;
    end
      plot(PixelsToVisualAngle([rSet; rSet; rSet]'), 100.*AccuracyOverEcc')
        legend('\beta=1.6', '\beta=1.65', '\beta=1.7')
        hold all;
        plot(PixelsToVisualAngle([rSet; rSet; rSet]'),100.*TA2Plot')
        xlabel('r'); ylabel('% correct');
    
end

% disp('Stop')
% for r=1:9
%     figure;
%     for beta = 1:3
%         clear tmp
%         tmp(:,:) = mean(Accuracy(:,r,:,:));
%         for p=1:4
%             phires(p) = mean([tmp(p,:) tmp(9-p,:)]);
%         end
%         plot(phires)
%     end
% end


figure('Position',[1 1 400 400])
% csvwrite(['acc' int2str(Pctr) '.csv'] , AccuracyOverEcc');
plot(PixelsToVisualAngle([rSet; rSet; rSet]'), 100.*mean(PersonAcc,3)')
legend('\beta=1.6', '\beta=1.65', '\beta=1.7')
hold all;
% plot(PixelsToVisualAngle([rSet; rSet; rSet]'),100.*TA2Plot', ':x')
xlabel('r'); ylabel('% correct');

linearmodel = -597+409.145*repmat([1.6,1.65,1.7],9,1)-11.441*repmat(PixelsToVisualAngle(rSet)',1,3);
linearmodel(linearmodel<0)=0;
plot(PixelsToVisualAngle([rSet; rSet; rSet]'), linearmodel)







function d = getDistances(Fixations)
x = (Fixations.x-512).^2;
y = (Fixations.y-512).^2;
d = sqrt(x+y);
end

end