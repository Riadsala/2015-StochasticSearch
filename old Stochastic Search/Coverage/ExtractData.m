function ExtractData
%%%% sort data for roughness experiment
clear all
close all


NumOfFixations = cell(5,3);
FinalFix2Target = cell(5,3);
NumOfFixationsTA = cell(5,3);


People = {'gw', 'hw', 'lm', 'jf', 'ps'};



parfor Pctr = 1:5
    Person = People{Pctr};
    [Results timeoffsets] = SortOutEventData(Person);
    for run = 1:3
        FixFile = fopen([Person int2str(run) 'FXD.txt']);
        data = fscanf(FixFile, '%c',inf);
        fclose(FixFile); 
        %%% delete header
        cursor = regexp(data, 'Fix number', 'once');
        data(1:(cursor-1)) = [];
        cursor = regexp(data, '1', 'once');
        data(1:(cursor-1)) = [];  
        % convert to matrix
        Fixs = str2num(data);
        for ctr = 1:100
            t = ctr+(run-1)*100;
            if Results(t).basicinfo(3)==1  
                % for each trial, get fix info
                StartTime = Results(t).basicinfo(1)-timeoffsets(run);
                EndTime = Results(t).basicinfo(2)-timeoffsets(run);
                Fixations = getFixInfo(StartTime, EndTime, Fixs);
                Results(t).basicinfo(9) = Fixations.number;
                Results(t).fixations = [Fixations.x, Fixations.y];   
            end
        end
    end
    P(Pctr).Results = Results;    
end
save 'Results.mat' P;   
end