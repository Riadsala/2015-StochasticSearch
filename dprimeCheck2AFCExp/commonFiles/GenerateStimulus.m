function  [stim x y clipping]  = GenerateStimulus(targPresent, a, tp1, tp2, useMask)
n = 1024; % size of search area
targetSpace = 'polar';
gaborSize = 50;
if nargin<2
    % dynamic range of noise
    a = 0.75; % default value;
end
if nargin < 4
    targetSpace = 'cart';
    tp1 = randi(1024-gaborSize);
    tp2 = randi(1024-gaborSize);
    
    useMask = false;
end

beta = 1; % roughness
alpha = 0.3; % strengh of target

im = GenerateFractals(n, beta,a);

if targPresent
    g = GenGabor(gaborSize);
    [stim x y] = PlaceTarget(im, g, alpha, targetSpace ,tp1, tp2);
else
    stim = im;
    x = NaN;
    y = NaN;
end

if useMask
    %% make mask
    X = repmat(1:n, [n,1])-n/2;
    Y = X';
    D =sqrt(X.^2 + Y.^2);
    D(D<n/2) = 1;
    D(D>1) = 0;
    stim = 255*(D .* stim+0.5);
else
    stim = 255*(stim+0.5);
end
clipping = sum(stim(:)<0) + sum(stim(:)>255);
end

function [out, x, y] = PlaceTarget(im, g, alpha, placement, loc1, loc2)
N = size(im,1);
if strcmp(placement, 'cart')
    
    x = loc1;
    y = loc2;
elseif strcmp(placement, 'polar')
    phi = loc1;
    r = loc2;
    x = round(r*cos(phi) + N/2-size(g,1)/2);
    y = round(r*sin(phi) + N/2-size(g,2)/2);
end

out = im;
out(x:(x+size(g,1)-1), y:(y+size(g,2)-1)) = out(x:(x+size(g,1)-1), y:(y+size(g,2)-1)) + alpha*g;
% correct [x,y] so we output the centre point, not the corner
x = x+size(g,1)/2;
y = y+size(g,2)/2;
end



function out = GenerateFractals(n, beta,a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function to generate a fractal surface of n x n
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%
% GENERATE SPECTRUM
%%%%%%%%%%%%%%%%%%%%%
% Generate indices
nq=n/2; V=-repmat((-nq:nq-1)', 1,2*nq); U=repmat((-nq:nq-1), 2*nq,1); f=sqrt(U.*U+V.*V);
theta=rand(n,n)*2*pi; % Generate random phase
mag=power(f,-beta); % Generate magnitude
mag=ifftshift(mag);  mag(1,1)=0; %shift & zero d.c.
[x,y]=pol2cart(theta,mag); F=x+1i*y; %convert to cartesian


out=ifft2(F, 'symmetric');
out = out - mean(out(:));
out = (0.5*a)* out/max(abs(out(:)));
% adjust so mean is 0

end




