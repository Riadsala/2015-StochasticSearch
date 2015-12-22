% function SimulatedSearch
clear all
close all
N = 1024;


numfix = 2^11

init_fix = [N/2, N/2];
c = [0 0];

SaccadeAmps = cumsum(csvread('SaccadeAmps.csv'));
SaccadeAmps(15:20) = [];
saccdists = [];

fix = zeros(numfix,2);
hotspot = zeros(1024);

fix(1,:) = init_fix;
fixhist = [fix];
for t=1:numfix
    hotspot(fix(t,1),fix(t,2)) = hotspot(fix(t,1),fix(t,2))+1;
    % make a random jump.  choose saccade amplitude
    p = rand*max(SaccadeAmps);
    saccampmin = VisualAngleToPixels(find(SaccadeAmps>p, 1 )/2-0.25);
    saccampmax = VisualAngleToPixels(find(SaccadeAmps>p, 1 )/2+0.25);
    X = repmat((1:N)-fix(t,1), [N,1]);
    Y = repmat((1:N)'-fix(t,2), [1,N]);
    dists = sqrt(X.^2+Y.^2);
    dists = (dists<saccampmax).*(dists>saccampmin);
    if sum(dists(:))==0
        disp('oops')
    end
    p = sum(dists(:))*rand;
    dists = cumsum(dists(:));
    [fix(t+1,1) fix(t+1,2)] = ind2sub(N,find(dists<p, 1, 'last' ));
    fixsaccamp(t) = norm(fix(t,:)-fix(t+1,:));
    ydiff(t) = fix(t,1)-fix(t+1,1);
end

% figure
% plot(fix(:,1),fix(:,2), 'r-x')
% axis([1 1024, 1 1024]);

hotspot = imfilter(hotspot, fspecial('Gaussian',13, 5));
figure;
imshow(hotspot,[]);


figure;
subplot(2,1,1)
fsa(1,:) = fft(ydiff(1:512));
fsa(2,:) = fft(ydiff(513:2*512));
fsa(3,:) = fft(ydiff(1025:3*512));
fsa(4,:) = fft(ydiff(1537:2048));
fsa(5,:) = fft(ydiff((1:512)+256));
fsa(6,:) = fft(ydiff((513:2*512)+256));
fsa(7,:) = fft(ydiff((1025:3*512)+256));

pfsa  = real(fsa).^2+imag(fsa).^2;
ps = mean(pfsa,1);
ps(1)=0;
plot(log(1:256), log(ps(1:256)))
subplot(2,1,2);

fsa(1,:) = fft(fixsaccamp(1:512));
fsa(2,:) = fft(fixsaccamp(513:2*512));
fsa(3,:) = fft(fixsaccamp(1025:3*512));
fsa(4,:) = fft(fixsaccamp(1537:2048));
fsa(5,:) = fft(fixsaccamp((1:512)+256));
fsa(6,:) = fft(fixsaccamp((513:2*512)+256));
fsa(7,:) = fft(fixsaccamp((1025:3*512)+256));

pfsa  = real(fsa).^2+imag(fsa).^2;
ps = mean(pfsa,1);
ps(1)=0;
plot(log(1:256), log(ps(1:256)))
hold on;
x=4.5:0.1:5.5;
y = 18+9-2*x;
plot(x,y,'r-')

x=3:0.1:4.5;

y = 19+2.5-x;
plot(x,y,'g-')
% end




