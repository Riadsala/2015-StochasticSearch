function getSaccadeDirection
close all
P = csvread('FourierComponents.csv');


theta = 0:0.01:2*pi;
size(theta)
cumv_q = P(1)*theta;
for j=2:64;
    cumv_q =  cumv_q + (-1)^(j-1).*real(P(j)).*sin((j-1)*theta)./(j-1) + (-1)^j*imag(P(j)).*cos((j-1)*theta)./(j-1);
end
cumv_q=cumv_q./64;

a = max(cumv_q);
phi = [];
for t=1:30000;
    x = a*rand;
    phi =[phi, theta(find(x<cumv_q, 1 ))];
end


figure;
theta=(0:pi/32:2*pi);
p2 = histc(phi,theta);
p2(65) = [];
theta(65) = [];
rose(phi)