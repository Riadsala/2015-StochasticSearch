close all


load searchSimResults

% read in target locations
targInfo = csvread('../analysis/trialInfo.txt', 1);


Beta = [1.60 1.65 1.70];

fout = fopen('searchSimResults.txt', 'w');
fprintf(fout, 'strat, beta,  nfix\n');

 fout2 = fopen('scanpaths.txt', 'w');
 fprintf(fout2, 'strategy, trial, beta, x,y, n \n');
t = 0;

            for b = 1:3
                for i = 1:length(res.ran.F)
                 t = t+1;
               
                beta = Beta(b);
                fprintf(fout, 'stochastic, %f, %d\n', beta, length(res.ran.F{b, i}));
                 fprintf(fout, 'optimal, %f, %d\n', beta, length(res.opt.F{b, i}));
                
                for f = 1:size(res.ran.F{b,i},1)
                    fprintf(fout2, 'stochastic, %d, %.3f, %d, %d, %d \n', ...
                        t, beta, f, res.ran.F{b,i}(f,1),res.ran.F{b,i}(f,2));
                end
                for f = 1:size(res.opt.F{b,i},1)
                    fprintf(fout2, 'optimal, %d, %.3f, %d, %d, %d \n', ...
                        t, beta, f, res.opt.F{b,i}(f,1),res.opt.F{b,i}(f,2));
                end
            end
        end

fclose(fout)



