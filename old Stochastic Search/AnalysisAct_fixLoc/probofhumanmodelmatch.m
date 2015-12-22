
% Area of circle containing all the model and the human saccade target for
% current fixation
r=2:0.5:8.5;
Area = pi*r.^2;
% best case area covered by model fixations
ma = 3*pi;
prob = ma./Area;
plot(r,prob, 'b-')
xlabel('radius of circle containing all possible saccade targets');
ylabel('probability of landing withing 1deg of a model point');
hold on;
% worst case... models potential saccade targets overlab (can't be this bad
% in practise).
ma = pi;
prob = ma./Area;
plot(r,prob, 'r-')
legend('Best case: models saccade targets do not overlap', 'worse than worst case - only one model saccade target')
plot(r,0.177,'k--')


