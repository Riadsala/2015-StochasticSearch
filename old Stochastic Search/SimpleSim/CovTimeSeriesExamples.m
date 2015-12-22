clear all
close all
load
figure;

for p=1:7
    subplot(2,4,p)
    clear idx
    idx = find(isfinite(CoverageTimeSeries(40,:,p)));
    ex = CoverageTimeSeries(:,idx(3),p);
    plot(ex);

end
      trialcoverage = csvread('randomwalktrialcoverage.csv');
      figure
      for p=1:8
         subplot(2,4,p);
         plot(trialcoverage(p+10,:), 'r');
      end
      
       trialcoverage = csvread('randomcoordstrialcoverage.csv');
      figure
      for p=1:8
         subplot(2,4,p);
         plot(trialcoverage(p+10,:), 'r');
      end
      