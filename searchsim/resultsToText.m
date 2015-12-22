close all
clear all

load searchSimResults

reps = 5;
R = [100 225 350];
Phi = 0:45:315;
Beta = [1.60 1.65 1.70];

fout = fopen('searchSimResults.txt', 'w');
fprintf(fout, 'strat, beta, ecc, phi, nfix\n');

 fout2 = fopen('scanpaths.txt', 'w');
 fprintf(fout2, 'strategy, trial, beta, ecc, phi, n, x, y \n');
t = 0;
for i = 1:reps
    for j = 1:length(R)
        for k = 1:length(Phi)
            for b = 1:3
                 t = t+1;
                [i j k b]
                r = R(j);
                phi = Phi(k);
                beta = Beta(b);
                fprintf(fout, 'stochastic, %f, %f, %f, %d\n', beta, r, phi, res.ran.nFix(i,j,k,b));
                fprintf(fout, 'optimal, %f, %f, %f, %d\n', beta, r, phi, res.opt.nFix(i,j,k,b));
                
                for f = 1:size(res.ran.F{i,j,k,b},1)
                    fprintf(fout2, 'stochastic, %d, %.3f, %d, %d, %d, %d, %d\n', ...
                        t, beta, r, phi, f, res.ran.F{i,j,k,b}(f,1),res.ran.F{i,j,k,b}(f,2));
                end
                for f = 1:size(res.opt.F{i,j,k,b},1)
                    fprintf(fout2, 'optimal, %d, %.3f, %d, %d, %d, %d, %d\n', ...
                        t, beta, r, phi, f, res.opt.F{i,j,k,b}(f,1),res.opt.F{i,j,k,b}(f,2));
                end
            end
        end
    end
end 
fclose(fout)



