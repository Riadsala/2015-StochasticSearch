
clear all
close all


N = 1024;
A = 0.05:0.002:0.15
a = 0;
for alpha = A
a = a+1;
for beta=1:3
 
    for r = 1:100
        fix = runSimRandom2(N, beta, [NaN NaN], alpha, 1, false);
      
       if isfinite(fix(2,1))
            fp(beta,a,r) = 1;

        else
            
            fp(beta,a,r) = 0;
        end
    end

end
end
plot(A ,100*mean(fp,3))
xlabel('alpha (threshold)');
ylabel('false positive rate');

export_fig modelFalsePosRates.pdf