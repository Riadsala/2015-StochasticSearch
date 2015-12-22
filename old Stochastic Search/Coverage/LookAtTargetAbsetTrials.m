clear all
load Results.mat


parfor Pctr=1:5
    R=P(Pctr).Results;
    ResultsTable = zeros(150,4);
    TAcounter = 0;
    for t=1:300
        % only look at TA trials
        if (R(t).basicinfo(3) == 1)
            TAcounter = TAcounter + 1;
            % check that trial was not false positive
            if (R(t).basicinfo(7) == 2)
                % trial was fasle positive, so report back NaNs
                ResultsTable(TAcounter,:) = [R(t).basicinfo(4),...
                    R(t).basicinfo(4), NaN, NaN];
            else
                if (R(t).basicinfo(9) > 0)
                    c = coverage(R(t).fixations, 1024);
                else
                    c=0;
                end
                ResultsTable(TAcounter,:) = [R(t).basicinfo(4),...
                    R(t).basicinfo(4), R(t).basicinfo(9), c];
            end
        end
    end
    Results(Pctr).Table = ResultsTable; %#ok<AGROW>
end
clear P
save Results2.mat