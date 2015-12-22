load ../Data/s7/Results


fout = fopen('parsedResults.txt', 'w');
fprintf(fout, 'beta phi ecc p n\n');


for b = 1:length(dat)    
  
    p = mean(dat(b).res);
   fprintf(fout, '%.3f %.3f %.3f %.3f %d\n', dat(b).beta, dat(b).phi, dat(b).ecc, p, 40);
end
fclose(fout);
