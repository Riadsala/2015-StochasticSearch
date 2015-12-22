function GenerateImageLibrary
tic
addpath('commonFiles');

% constants for outputting rendered image
slant = 60; tilt = 90;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta_Set = [1.6, 1.65, 1.7];
RMS = 1.1;
N = 1024;
NumberOfExamplesPerComb = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defect Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
targa = 10; targb = 10; c = 2;
B= [0 1 0; 1 2 1; 0 1 0];
vol_threshold = 50;


nAngles = 8;
nEccs = 4;
nBlocks = 2;

% make conditions
ctr = 0;
for blk = 1:nBlocks
    for b = 1:length(beta_Set);
        for phi = 1:nAngles
            for ecc = 1:nEccs
                ctr = ctr + 1;
                dat(ctr).beta = beta_Set(b);
                dat(ctr).phi = (phi/nAngles) * (2*pi);
                dat(ctr).ecc = ecc*(N)/(2*(nEccs+1));
                dat(ctr).res = [];
            end
        end
        ctr = ctr  +1;
        dat(ctr).beta = beta_Set(b);
        dat(ctr).phi = 0;
        dat(ctr).ecc = 0;
        dat(ctr).res = [];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GenerateFractals(N, beta_Set, NumberOfExamplesPerComb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add defects to textures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(dat)
    for ctr = 1:NumberOfExamplesPerComb
         r= dat(i).ecc;
         phi = dat(i).phi;
                T = csvread([int2str(ctr) '.csv']);
                T = T.*RMS;                             
                
                filename = ['ims/tp_beta=' num2str(dat(i).beta) '_r=' int2str(r) '_phi=' int2str(180*phi/pi) '_' int2str(ctr)];
                % generate random coordinates for defect placement
                [x,y] = pol2cart(phi, r)
                 x = round(x+N/2);
                y = round(y+N/2);
                [phi, r]
               
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % create defect
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                D = GenerateEllipsoidDefect(targa,targb,c);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % cut defect out of surface
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                T = AddDefect2Surface(T, D,vol_threshold,x,y,c);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % render and apply self shadowing
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                R = lambertian(T,slant,tilt, N);
                R = (R>0).*R;
                imwrite(R, [filename '.png']);
           
    end
end
toc
end




function GenerateFractals(n, beta_Set, NumberOfExamplesPerComb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function to generate a fractal surface of n x n
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slant = 60; tilt = 90;
N=1024
seed = 0;
for beta = beta_Set
    for ctr = 1:NumberOfExamplesPerComb
        seed = seed+1;
        %%%%%%%%%%%%%%%%%%%%%
        % GENERATE SPECTRUM
        %%%%%%%%%%%%%%%%%%%%%
        % Generate indices
        nq=n/2; V=-repmat((-nq:nq-1)', 1,2*nq); U=repmat((-nq:nq-1), 2*nq,1); f=sqrt(U.*U+V.*V);
        rand('state',seed); theta=rand(n,n)*2*pi; % Generate random phase
        mag=power(f,-beta); % Generate magnitude
        mag=ifftshift(mag);  mag(1,1)=0; %shift & zero d.c.
        [x,y]=pol2cart(theta,mag); F=x+1i*y; %convert to cartesian
        
        ht=ifft2(F, 'symmetric');
        % adjust so mean is 0
        ht = ht - mean(ht(:));
        ht = ht./RMSroughness(ht);
        basename=[int2str(seed) '.csv'];
        csvwrite(basename, ht);
        
        % save TA version
          R = lambertian(ht,slant,tilt, N);
          R = (R>0).*R;
          imwrite(R, ['ims/ta_beta' num2str(beta) '_' int2str(ctr) '.png']);
        
    end
end

end





function D = GenerateEllipsoidDefect(a,b,c)
a2 = a.*a; b2 = b.*b; c2 = c.*c;
for row=-2*a:2*a
    row_e=row+2*a+1;
    for col=-2*b:2*b
        col_e=col+2*b+1;
        x_dash=row;
        y_dash=col;
        sqr=c2*(1-(x_dash*x_dash/a2+y_dash*y_dash/b2));
        if sqr>0
            D(row_e, col_e)=sqrt(sqr);
        else
            D(row_e, col_e)=0;
        end
    end
end
B= [0 1 0; 1 2 1; 0 1 0];
D = filter2(B, D)./6;
end


function Td = AddDefect2Surface(T, D, vol_threshold, x, y, c)
step = 0.1;
n = size(D);
Td = T;
localmaxheight = max(max(T(x:x+n(1), y:y+n(2))));


height = c+localmaxheight;
D = height - D;
vol = 0;
while vol<vol_threshold
    % add defect to texture
    for i = 1:n(1)
        for j = 1:n(2)
            if D(i,j) ~= D(1,1)
                Td(i+x-floor(n(1)/2), j+y-floor(n(2)/2)) = min(Td(i+x-floor(n(1)/2), j+y-floor(n(2)/2)), D(i,j));
            end
        end
    end
    % check difference
    vol = sum(sum(T - Td));
    D = D - step;
end
end