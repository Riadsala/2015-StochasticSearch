close all
SA = csvread('SaccadeAngles.csv');
n=64;
SA=mod(SA,2*pi);
theta=(0:(2*pi/n):2*pi);

p = histc(SA,theta);
p(n+1)=[]
theta(n+1)=[]
bar(theta, p./a, 'histc')
% figure

% subplot(2,1,1);
a=trapz(theta,p);
% plot(theta, p/a, 'b');

hold on;


P = fft(p/a);


q = P(1);
for j=2:n;
    q =  q + (-1)^(j-1).*real(P(j)).*cos((j-1)*theta) + (-1)^(j-1)*imag(P(j)).*sin((j-1)*theta);

end
q=q./n;


close all
figure
p(65)=p(1);
q(65) = q(1);
theta(65) = 2*pi;
polar(theta,p/a, 'b')
hold on
polar(theta,q, 'r')
legend('Empirical Data', 'Fourier Component Data', 'location', 'Southwest');
%  q = circshift(q, [1 -1]);


theta(65) = 2*pi;

cumv_q = 1/n*(P(1)*theta...
    - real(P(2)).*sin(theta) + imag(P(2)).*cos(theta)...
    + real(P(3)).*sin(2*theta)/2 - imag(P(3)).*cos(2*theta)/2 ...
    - real(P(4)).*sin(3*theta)/3 + imag(P(4)).*cos(3*theta)/3 ...
    + real(P(5)).*sin(4*theta)/4 - imag(P(5)).*cos(4*theta)/4 ...
    - real(P(6)).*sin(5*theta)/5 + imag(P(6)).*cos(5*theta)/5 ...
    + real(P(7)).*sin(6*theta)/6 - imag(P(7)).*cos(6*theta)/6 ...
    - real(P(32)).*sin(31*theta)/31+ imag(P(32)).*cos(31*theta)/31 ...
    + real(P(31)).*sin(30*theta)/30 - imag(P(31)).*cos(30*theta)/30 ...
    - real(P(30)).*sin(29*theta)/29 + imag(P(30)).*cos(29*theta)/29 ...
    + real(P(29)).*sin(28*theta)/28 - imag(P(29)).*cos(28*theta)/28 ...
    - real(P(28)).*sin(27*theta)/27 + imag(P(28)).*cos(27*theta)/27 ...
    + real(P(27)).*sin(26*theta)/26 - imag(P(27)).*cos(26*theta)/26);



% a=max(cumv_q);
% subplot(2,1,2);
% plot(abs(P))

csvwrite('FourierComponents.csv', P);
cumv_q(1)=0;
% simulate random sampling.
phi = [];
for t=1:30000;
    x = a*rand;
    phi =[phi, theta(find(x<cumv_q, 1 ))-pi/64*rand];
end


% figure;
% theta=(0:pi/32:2*pi);
% 
% p2 = histc(phi,theta);
% p2(65) = [];
% theta(65) = [];
% plot(theta,p/sum(p), 'b')
% hold on
% plot(theta,p2/sum(p2), 'r')


