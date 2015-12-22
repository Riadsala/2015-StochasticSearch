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

Person = {'LM', 'MY'};
preceedTrialStartTime = 50;
postTrialEndTime = -100;
maxAllowedDist=60;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sort out EventData - filenames and keypresses
%%%%%%help %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for Pctr=1:2
%     for part = 1:20
%        
%         FixFile = fopen( [Person{Pctr} int2str(part) 'EFD.txt']);
% 
% 
%         FixData = fscanf(FixFile, '%c',inf);
%         fclose(FixFile);
%         %%% delete header
%         cursor = regexp(FixData, 'Timestamp', 'once');
%         FixData(1:(cursor-1)) = [];
%         cursor = regexp(FixData, '\d', 'once');
%         FixData(1:(cursor-1)) = [];
%         % convert to matrix - first need to replace validity strings with
%         % numbers
%         FixData=regexprep(FixData, 'Both', '1');
%         FixData=regexprep(FixData, 'None', '0');
%         FixData=regexprep(FixData, 'LeftOnly', '0');
%         FixData=regexprep(FixData, 'RightOnly', '0');
%         Fixs = str2num(FixData);
%         clear FixData
%         EventFile = fopen([Person{Pctr} int2str(part) 'EVD.txt']);
%         EventData = fscanf(EventFile, '%c',inf);
%         fclose(EventFile);
%         %%% delete header
%         cursor = regexp(EventData, 'Description', 'once');
%         EventData(1:(cursor+14))=[];
%         N = size(EventData,2);
%         % find fixcross.jpg.
%         FixcrossLocations = regexp(EventData, 'fixcross.jpg');
%         %     n = size(FixcrossLocations, 2);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % go through all the trials in the run and get break times
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         ImageLocations = regexp(EventData, '_beta');
%         T = size(ImageLocations,2);
%         for b = 1:BlocksInRun
%             for t = 1:TrialsInBlock
%                 TrialCounter = TrialCounter+1;
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % find relvant fixation cross and get following image display time
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 cursor = FixcrossLocations((b-1).*TrialsInBlock+t);
%                 TrialStartTime = str2double(regexp(EventData(cursor:N), '\s\d+\s', 'match', 'once'))-preceedTrialStartTime;
%                 cursor = ImageLocations((b-1).*TrialsInBlock+t);
%                 TrialEndTime = str2double(regexp(EventData(cursor:N), '\s\d+\s', 'match', 'once'))+postTrialEndTime;
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % Check fixation
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Fidx = find((Fixs(:,1)<TrialEndTime).*(Fixs(:,1)>TrialStartTime));
%                 Fixations.y = 1024-(1200-Fixs(Fidx, 4)-xoffset);
%                 Fixations.x = Fixs(Fidx, 3)-yoffset;
%                 DistancesFromCentre = getDistances(Fixations);
%                 Results(TrialCounter, 5) = mean(DistancesFromCentre);
%                 Results(TrialCounter, 6) = std(Fixations.x);
%                 Results(TrialCounter, 7) = std(Fixations.y);
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % get surface parameters from filename
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 cursor = ImageLocations((b-1).*TrialsInBlock+t);
%                 Filename = regexp(EventData(cursor:N), 'beta[._=\w]*png',  'once', 'match');
%                 beta = str2double(regexp(Filename, '(?<=beta=)1\.[567]*', 'match'));
%                 r = str2double(regexp(Filename, '(?<=r=)(50|100|150|200|250|300|350|400|450)', 'match'));
%                 if r>0
%                     phi = str2double(regexp(Filename, '(?<=phi=)(0|45|90|135|180|225|270|315)', 'match'));
%                 else
%                     r = inf;
%                     phi = inf;
%                 end
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % new filter - include trials with mean,std fix < threshold
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 dRule = (Results(TrialCounter, 5)<60)*(Results(TrialCounter, 6)<40)*(Results(TrialCounter, 7)<40);
%                 if dRule&(size(Fidx,1)>0)
%                     % get Right|Left Keyboard response
%                     Response =  regexp(EventData(cursor:N), '(Right)|(Left)',  'once', 'match');
%                     % check that response was found before the next trial was shown
%                     % if not, set response to 2 (do not count)
%                     if (b-1)*TrialsInBlock+t~=T
%                         RespLoc =  regexp(EventData(cursor:N), '(Right)|(Left)',  'once');
%                         if (RespLoc+cursor)>ImageLocations((b-1).*TrialsInBlock+t+1)
%                             %                         disp('problem with missing response');
%                             %                         disp(part)
%                             %                         disp(Filename)
%                             Results(TrialCounter,4) = 2;
%                         else
%                             if strcmp(Response, 'Right')
%                                 Results(TrialCounter,4) = 0;
%                             elseif strcmp(Response, 'Left')
%                                 Results(TrialCounter,4) = 1;
%                             end
%                         end
%                     end
%                 elseif size(Fidx,1)==0
%                     Results(TrialCounter,4) = 3;
%                 else
%                     Results(TrialCounter,4) = 2;
%                 end
%                 Results(TrialCounter,1:3) = [beta, r, phi];
%             end
%         end
%     end
% end
% csvwrite('results.csv', Results);
Results = csvread('results.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sort out Results - how does roughness and r effect signal detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
betaSet = [1.6, 1.65, 1.7];
rSet = 50:50:450;
rSet = rSet;
phi = [0,45,90,135,180,225,270,315];
Accuracy = zeros(3,8,8);
AccuracyOverEcc = zeros(3,9);
TargetAbsentAccuracy = zeros(3,1);
TA2Plot = zeros(3,9);
for b = 1:3
    for r=1:9
        for p=1:8
            beta = betaSet(b);
            index = find((Results(:,1)==beta).*(Results(:,2)==rSet(r)).*(Results(:,3)==phi(p)).*(Results(:,4)~=2));
            Accuracy(b,r,p) = sum(Results(index, 4)==1)./size(Results(index, 4)==1, 1);
        end
        tmp = Accuracy(b,r,:);
        tmp=tmp(:);
        AccuracyOverEcc(b,r) = mean(tmp(isfinite(tmp)));
    end
    index = find((Results(:,1)==beta).*isinf(Results(:,2)).*(Results(:,4)~=2));
    TargetAbsentAccuracy(b) = 1-mean(Results(index, 4)); %#ok<FNDSB>
    TA2Plot(b,:) = TargetAbsentAccuracy(b)*ones(1,9);
end
AccuracyOverEcc;
csvwrite('LMacc.csv', AccuracyOverEcc');
plot(PixelsToVisualAngle([rSet; rSet; rSet]'), 100.*AccuracyOverEcc', '-x')
legend('\beta=1.6', '\beta=1.65', '\beta=1.7')
hold all;
plot(PixelsToVisualAngle([rSet; rSet; rSet]'),100.*TA2Plot', ':x')
xlabel('r'); ylabel('% correct');

dprime = csvread('dprime.csv');
figure;
% dprime =  (icdf('norm', AccuracyOverEcc,0, 1)-icdf('norm',1- TA2Plot,0, 1))';
plot(PixelsToVisualAngle([rSet; rSet; rSet]'),dprime)
title('Effect of eccentricity and \beta on d'' ');
xlabel('r  degrees')
ylabel('d''');
csvwrite('dprime.csv', dprime);


% simple multilinear regression
figure;
plot(PixelsToVisualAngle([rSet; rSet; rSet]'), 100*AccuracyOverEcc', '-x')
hold on
 linearmodel = -596.5+409.145*repmat([1.6,1.65,1.7],9,1)-11.441*repmat(PixelsToVisualAngle(rSet)',1,3);
linearmodel(linearmodel<0)=0;
plot(PixelsToVisualAngle([rSet; rSet; rSet]'), linearmodel, ':x')
xlabel('r (visual angle)'); ylabel('% correct');
% title('multilinear regression. R squared=0.942')

figure;
plot(PixelsToVisualAngle([rSet; rSet; rSet]'), dprime, '-x')
hold on
linearmodel = repmat([1.387, 2.877, 3.941],9,1)-repmat([0.004, 0.006, 0.007],9,1).*repmat(rSet',1,3);

plot(PixelsToVisualAngle([rSet; rSet; rSet]'), linearmodel, ':x')
title('linear regression in just r - do each \beta seperatly. Rsq = [0.728, 0.981, 0.972]')








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get number of scrapped trials
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PropOfTrialsScrapped = 100*sum(Results(:,4)>=2)./size(Results,1);
disp([int2str(PropOfTrialsScrapped) '% of all trials have been scrapped'])

 model = -4.205+repmat([1.6, 1.65, 1.7]',1,9).*3.056+repmat(50:50:450,3,1)*(-0.002);
 model(model<0)=0;
 plot(50:50:450, 100*model, '-')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % do something with phi data - average over r
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accuracy = zeros(3,8,9);
% AccuracyOverPhi = zeros(3,8);
% TargetAbsentAccuracy = zeros(3,1);
% TA2Plot = zeros(3,8);
% clear tmp
% figure
% hold on;
% for b = 1:3
%     for p=1:8
%         for r=1:9
%             beta = betaSet(b);
%             index = find((Results(:,1)==beta).*(Results(:,2)==rSet(r)).*(Results(:,3)==phi(p)).*(Results(:,4)~=2));
%             Accuracy(b,p,r) = sum(Results(index, 4)==1)./size(Results(index, 4)==1, 1);
%         end
%         tmp = Accuracy(b,p,:);
%         tmp=tmp(:);
%         AccuracyOverPhi(b,p) = mean(tmp(isfinite(tmp)));
%     end
%     index = find((Results(:,1)==beta).*isinf(Results(:,2)).*(Results(:,4)~=2));
%     TargetAbsentAccuracy(b) = 1-mean(Results(index, 4)); %#ok<FNDSB>
%     TA2Plot(b,:) = TargetAbsentAccuracy(b)*ones(1,8);
% 
%     AccuracyOverPhi(b,9) = AccuracyOverPhi(b,1);
%     [x y] =pol2cart(pi*([phi 360]-90)/180, AccuracyOverPhi(b,:));
%     plot(x,y, '-x')
% end
% 
% plot(0,0,'rx')
% plot([-0.8, 0.8], [0 0], '-k');
% plot([0 0], [-0.8, 0.8], '-k');
% axis([-0.8 0.8 -0.8 0.8])
% axis equal
% 
% Accuracy(:,9,:)= Accuracy(:,1,:);
% 
% for b=1:3
%     figure
%     hold all
%     for r=1:9
%         [x y] =pol2cart(pi*([phi 360]-90)/180, Accuracy(b,:,r));
%         plot(x,y, '-x')
%     end
%     plot(0,0,'rx')
%     plot([-1, 1], [0 0], '-k');
%     plot([0 0], [-1, 1], '-k');
%     axis([-1 1 -1 1])
%     axis equal
%     % legend('r=50','r=100','r=150','r=200','r=250','r=300','r=350','r=400','r=450')
% end


    function d = getDistances(Fixations)
        x = (Fixations.x-512).^2;
        y = (Fixations.y-512).^2;
        d = sqrt(x+y);
    end

end