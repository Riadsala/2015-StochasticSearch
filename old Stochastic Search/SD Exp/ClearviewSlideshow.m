function ClearviewSlideshow
close all
clear all

seed = 0;
overallctr = 0;
rand('state', seed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta_Set = [1.6, 1.65, 1.7];
r_Set = [50, 100, 150, 200, 250, 300, 350, 400, 450];
N = 1024;
NumberOfExamplesPerComb = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add defects to textures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defect_positions = zeros(N);
pctr=0;
for beta = beta_Set
    for ctr = 1:NumberOfExamplesPerComb
        seed = seed+1;
        for phi = 0:45:315
            for r = r_Set
                pctr = pctr + 1;
                filenames{pctr, :}  = ['_beta=' num2str(beta) '_r=' int2str(r) '_phi=' int2str(phi) '_seed=' int2str(seed) '.png'];
                params(pctr,:) = [beta, r, seed];
            end
            pctr = pctr+1;
            filenames{pctr,:} = ['_beta=' num2str(beta) '_r=TA' '_seed=' int2str(seed) '.png'];
            params(pctr,:) = [beta, NaN, seed];
            pctr = pctr+1;
            filenames{pctr,:} = ['_beta=' num2str(beta) '_r=TA' '_seed=' int2str(seed) '.png'];
            params(pctr,:) = [beta, NaN, seed];

        end
    end
end
numberoftrials = pctr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set seed for random slideshow order!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seed = 2;
rand('state', seed);
randn('state', seed);
ranindex = randperm(numberoftrials);
filenames = filenames(ranindex);
params = params(ranindex,:);
SaveData(filenames, numberoftrials);
end

function SaveData(filenames, numberoftrials)
tctr=0;
for run = 1:20
    fid      = fopen(['SlideshowIndex' int2str(run) '.txt'],'w');
    for batch = 1:((numberoftrials/20)/4)
        fprintf(fid,'%s\t', 'breakrest.jpg');
        fprintf(fid,'%s\t', 'keypress'); % how many seconds to display the fixation cross
        fprintf(fid,'%s\t', 'white');    % background colour
        fprintf(fid,'\r\n');
        for trial = 1:4
            tctr=tctr+1;
            fprintf(fid,'%s\t', 'fixcross.jpg');
            fprintf(fid,'%s\t', '2'); % how many seconds to display the fixation cross
            fprintf(fid,'%s\t', 'white');    % background colour
            fprintf(fid,'\r\n');
            fprintf(fid,'%s\t', strcat(filenames{tctr})); %stimulus filename
            fprintf(fid,'%s\t', '0.2'); %tell Clearview to wail til keypress, otherwise specify how many seconds to wait for
            fprintf(fid,'%s\t', 'white');
            fprintf(fid,'\r\n');
            fprintf(fid,'%s\t', 'mask.jpg');
            fprintf(fid,'%s\t', '0.5'); % how many seconds to display the fixation cross
            fprintf(fid,'%s\t', 'white');    % background colour
            fprintf(fid,'\r\n');  
            fprintf(fid,'%s\t', 'respdot.jpg');
            fprintf(fid,'%s\t', 'keypress'); % how many seconds to display the fixation cross
            fprintf(fid,'%s\t', 'white');    % background colour
            fprintf(fid,'\r\n');
        end
        
    end
    fclose(fid);
end
end