function SimulatedSearch
clear all
close all
N = 1024;

init_fix = [N/2, N/2];
c = [0 0];

SaccadeAmps = cumsum(csvread('SaccadeAmps.csv'));
SaccadeAmps(15:20) = [];
saccdists = [];
betaS = [1.6, 1.65, 1.7];
rvalues = [100 225 350];
for beta = 1:3
    beta
    parfor rcase = 1:3
        results = [];
        r=rvalues(rcase);
        phi = rand*360;
        [x y] = pol2cart(phi, r);
        x = ceil(x); y = ceil(y);
        targetloc = [x y]+N/2;
        for trial=1:100
            targetfound = 0;
            numfix = 0;
            fix = init_fix;
            fixhist = [fix];
            while targetfound==0
                oldfix=fix;
%                 if numfix>200
%                     numfix = NaN;
%                     break
%                 end
                numfix=numfix+1;
                % check to see if target is detected, and if so, make saccade to target
                dist2target  = PixelsToVisualAngle(norm(fix-targetloc));
                probofdetection = linear_model(betaS(beta), dist2target, c);
                chance = rand*100;
                if chance<probofdetection
                    targetfound = 1;
                    fix = targetloc;
                else
                    % make a random jump.  choose saccade amplitude
                    p = rand*max(SaccadeAmps);
                    saccampmin = VisualAngleToPixels(min(find(SaccadeAmps>p))/2-0.25);
                    saccampmax = VisualAngleToPixels(min(find(SaccadeAmps>p))/2+0.25);
                    X = repmat((1:N)-fix(1), [N,1]);
                    Y = repmat((1:N)'-fix(2), [1,N]);
                    dists = sqrt(X.^2+Y.^2);
                    dists = (dists<saccampmax).*(dists>saccampmin);
                    if sum(dists(:))==0
                       disp('oops') 
                    end
                  
                    
                    p = sum(dists(:))*rand;
                    dists = cumsum(dists(:));
                    [fix(1) fix(2)] = ind2sub(N,find(dists<p, 1, 'last' ));

                end
                fixhist = [fixhist; fix];
                saccdists = [saccdists, norm(oldfix-fix)];
            end
            results(trial) = numfix;   
            %         figure
            %         plot(fixhist(:,1),fixhist(:,2), 'r-x')
            %         axis([1 1024, 1 1024]);

        end
        meansaccades(beta,rcase) = mean(results(isfinite(results)));
        stdsaccades(beta, rcase) = std(results(isfinite(results)));
    end
end
figure;
hist(PixelsToVisualAngle(saccdists),0.5:0.5:20);
csvwrite('meansacc.csv', meansaccades);

csvwrite('stdnsacc.csv', stdsaccades);  figure;
    errorbar(repmat([1.6,1.65,1.7]',1,3), meansaccades, stdsaccades./50);
    legend('r=150', 'r-225', 'r=350', 'location', 'northeast')

end

    function p = linear_model(beta, r, c)


        p = -597 + 409.145*beta - 11.441*r;

    end
