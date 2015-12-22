% function SimulatedSearch
clear all
close all
N = 1024;


numfix = 1000;

init_fix = [N/2, N/2];
c = [0 0];

weibull_a = 4.62;
weibull_b = 1.50;


fix = zeros(numfix,2);
hotspot = zeros(1024);

fix(1,:) = init_fix;
fixhist = [fix];
for t=1:numfix
    t
    hotspot(fix(t,1),fix(t,2)) = hotspot(fix(t,1),fix(t,2))+1;
    % calculate polar coords for all pixels
    %     X = repmat((1:N)-fix(t,1), [N,1]);
    %     Y = repmat((1:N)'-fix(t,2), [1,N]);
    %     dists = PixelsToVisualAngle(sqrt(X.^2+Y.^2));
    %     dists = wblpdf(dists,weibull_a, weibull_b);
    validfix=0;
    while validfix == 0
        dist = VisualAngleToPixels(wblinv(rand, weibull_a, weibull_b));
        direction = rand*2*pi;
        [x y] = pol2cart(direction, dist);
        x=round(x)+fix(t,1);
        y=round(y)+fix(t,2);
        validfix = (0<x).*(x<=1024).*(0<y).*(y<=1024);
    end
    fix(t+1,:) = [x y];
    fixsaccamp(t) = norm(fix(t,:)-fix(t+1,:));
    ydiff(t) = fix(t,1)-fix(t+1,1);
end

% hist(PixelsToVisualAngle(fixsaccamp), 0.5:1:20)
% figure
% plot(fix(:,1),fix(:,2), 'r-x')
% axis([1 1024, 1 1024]);
% figure
% hotspot = imfilter(hotspot, fspecial('Gaussian',13, 5));
% figure;
% imshow(hotspot,[]);


% figure;
% subplot(2,1,1)
% fsa(1,:) = fft(ydiff(1:512));
% fsa(2,:) = fft(ydiff(513:2*512));
% fsa(3,:) = fft(ydiff(1025:3*512));
% fsa(4,:) = fft(ydiff(1537:2048));
% fsa(5,:) = fft(ydiff((1:512)+256));
% fsa(6,:) = fft(ydiff((513:2*512)+256));
% fsa(7,:) = fft(ydiff((1025:3*512)+256));
%
% pfsa  = real(fsa).^2+imag(fsa).^2;
% ps = mean(pfsa,1);
% ps(1)=0;
% plot(log(1:256), log(ps(1:256)))
% subplot(2,1,2);
%
% fsa(1,:) = fft(fixsaccamp(1:512));
% fsa(2,:) = fft(fixsaccamp(513:2*512));
% fsa(3,:) = fft(fixsaccamp(1025:3*512));
% fsa(4,:) = fft(fixsaccamp(1537:2048));
% fsa(5,:) = fft(fixsaccamp((1:512)+256));
% fsa(6,:) = fft(fixsaccamp((513:2*512)+256));
% fsa(7,:) = fft(fixsaccamp((1025:3*512)+256));
%
% pfsa  = real(fsa).^2+imag(fsa).^2;
% ps = mean(pfsa,1);
% ps(1)=0;
% plot(log(1:256), log(ps(1:256)))
% hold on;
% x=4.5:0.1:5.5;
% y = 18+9-2*x;
% plot(x,y,'r-')
%
% x=3:0.1:4.5;
%
% y = 19+2.5-x;
% plot(x,y,'g-')
% % end
%
%
%
%
