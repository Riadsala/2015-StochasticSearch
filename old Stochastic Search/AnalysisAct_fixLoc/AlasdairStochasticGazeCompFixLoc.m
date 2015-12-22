function numOfFixs = AlasdairStochasticGaze(I,SalMap, params, Location, seed1, seed2, Fix)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alasdair's Saliency Alg
% Version Three - now checks itself if it has found the target
% Have also tidied up a lot.
%
% Version Two - Just do Gabor Filtering and max/median weighting
%
% What have I done?
% Removed contrast channel. The contrast channel is essentially a bank of
% bandpass filters. All the gabor channels for any given scale approximate
% a bandpass filter, so it loos like there is possible redundancy
% No point taking centre surrounds of orientation channel as that's jsut
% essentially applying a bandpass, and the gabor filters are already
% kind of like a bandpass filter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rand('state', seed1*seed2);
N=max(size(SalMap));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement Inhibition of Return (IOR) and Eccentricity Dependant
% Processing (EDP).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IOR_locations = zeros(params.IOR_decaytime,2); % this is used to store the locations where IOR is in effect.
IOR_mask = createIOR_mask;
% set initial fixation to image centre
currentFixation = [N./2 N./2];
TargetFixated = 0; numOfFixs = 0;
inputfilename = ['tmp'];
% now carry out search for target...
wrj = 0;
while (TargetFixated==0)
    numOfFixs = numOfFixs + 1;
    %                create inhibition of return map for this fixation
    IOR_map = createIOR_map;
    %                 apply exponential EDP mask and IOR mask
    tSalMap = SalMap.*EDPmask(currentFixation).*IOR_map;
    imwrite(SimplyNormalise(tSalMap), ['tSalMap' int2str(numOfFixs) '.png']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate new fixation - look at the 5 largest peaks in the saliency
    % map and choose one at random
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fSalMap = tSalMap;
    for f=1:3
        fixs(f,1:2) = max2(fSalMap);
        fixs(f,3) = tSalMap(fixs(f,1), fixs(f,2));
        fSalMap = fSalMap.*(1-circshift(IOR_mask, fixs(f,1:2)));
    end
    LocChoice = sum(fixs(:,3))*rand;   
    if LocChoice  <= fixs(1,3)
        currentFixation = fixs(1,1:2);
    elseif LocChoice <= (fixs(1,3)+fixs(2,3))
        currentFixation = fixs(2,1:2);
        wrj = wrj+1;
    elseif LocChoice <= (fixs(1,3)+fixs(2,3)+fixs(3,3))
        wrj = wrj+1;
        currentFixation = fixs(3,1:2);
    elseif LocChoice >  (fixs(1,3)+fixs(2,3)+fixs(3,3))
        disp('this should not happen');
    end
%     DisplayOutputMap(tSalMap, currentFixation, [inputfilename 'Fixation' int2str(numOfFixs) '.jpg']);
    %                 check if new fixations TargetFixatedes TargetLocation
    %                 [numOfFixs Fixation Location]
    TargetFixated = CheckForTargetFixated(currentFixation, Location);
    if numOfFixs > 500
        numOfFixs = NaN;
        break
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested Functions Live Below Here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function to create the IOR_mask
    function IOR_mask = createIOR_mask
        % use a gaussian distribution as a mask.
        IOR_mask = fspecial('gaussian',N, params.IOR_mask_size);
        % move gaussian to centre of image and normalise to [0, 1].
        IOR_mask = SimplyNormalise(circshift(IOR_mask, [-N/2 -N/2]));
    end
% function used to create the composite IOR map
    function IOR_map = createIOR_map
        % move all previous fixation locations back one step in time
        IOR_locations = circshift(IOR_locations, 1);
        % add current Fixation location to IOR list.
        IOR_locations(1, :) = currentFixation;
        %compute composite IOR map
        IOR_map = ones(N);
        for i = 1:min(numOfFixs, params.IOR_decaytime)
            IOR_map =IOR_map.*(1-(2/(i+1)).*circshift(IOR_mask, IOR_locations(i,:)));
        end
    end
% function used to create eccentricity dependant processing map.
    function EDP = EDPmask(Fixation)
        X = repmat((1:N)-Fixation(2), N,1);
        Y = repmat((1:N)'-Fixation(1),1,N);
        dist = sqrt(X.^2+Y.^2);
        EDP = exp(params.EDPradius.*dist);
    end
% function to check if the model is within params.FixationRadius of
% target's location
    function TargetFixated = CheckForTargetFixated(Fixation, Location)
        TargetFixatedx = (Fixation(1) > Location(1)-params.FixationRadius)&&...
            (Fixation(1) < Location(1)+params.FixationRadius);
        TargetFixatedy = (Fixation(2) > Location(2)-params.FixationRadius)&&...
            (Fixation(2) < Location(2)+params.FixationRadius);
        TargetFixated = TargetFixatedx * TargetFixatedy;
    end
end % end of main function


function DisplayOutputMap(map, loc, filename)
N = size(map,1);
map = SimplyNormalise(map);
loc1f = max(loc(1)-5,1);
loc1l = min(loc(1)+5,N);
loc2f = max(loc(2)-5,1);
loc2l = min(loc(2)+5,N);
output(:,:,1 ) = map;
output(:,:,2 ) = map;
output(:,:,3 ) = map;
output(loc1f:loc1l, loc2f:loc2l, 1)=output(loc1f:loc1l, loc2f:loc2l, 1)+0.5;
output(loc1f:loc1l, loc2f:loc2l, 2)=output(loc1f:loc1l, loc2f:loc2l, 1)-0.5;
output(loc1f:loc1l, loc2f:loc2l, 3)=output(loc1f:loc1l, loc2f:loc2l, 1)-0.5;
imwrite(output, filename);
end




function OriMap = GenerateOriMaps(I, params)
N=max(size(I));         %Get image dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first linear filtering stage - only use gabors gabors. Carry out
% filtering and normalisation and crossscale sum all in the same loop to
% save memory!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = fspecial('gaussian', N, N/16);
Ifouriertransform = fftshift(fft2(I));
% generate u and v indices - needed for creating filter
tmp = -N/2:N/2-1;
u = repmat(tmp, N, 1); u2 = u.^2;
v = repmat(tmp', 1, N); v2 = v.^2;
clear tmp
% create empty matrices to put filters and responses in later.
OriMap = zeros(N,N,params.NmbrOri);
% GaborFilters = zeros(1024, 1024);
for ori = 1:params.NmbrOri
    for lvl=1:params.NumberOfRes
        % generate filter in fourier domain
        % generate filter mask here rather than use Al_Gabor. i've managed
        % to speed things up a little by doing it this way.
        r = params.Gabor_r^(lvl+2+params.start_l);
        sigma2 = (params.Gabor_stddev^(lvl-1+params.start_l)).^2;
        phi = params.orientations(ori);
        u0 = r.*cosd(phi);   v0 = r.*sind(phi);
        bit1 = (u2 - 2*u0.*u + u0^2)/sigma2;
        bit2 = (v2 - 2*v0.*v + v0^2)/sigma2;
        G1= exp(-1/2*(bit1+bit2));
        bit1 = (u2 + 2*u0.*u + u0^2)/sigma2;
        bit2 = (v2 + 2*v0.*v + v0^2)/sigma2;
        G2= exp(-1/2*(bit1+bit2));
        GaborFilter = G1+G2;
%         GaborFilters = GaborFilters + GaborFilter;
        clear G1 G2 bit1 bit2
        % apply filter/mask to input image in fourier domain
        GaborResp =  real(ifft2(ifftshift(Ifouriertransform.*GaborFilter)));
        % squarewave rectify and normalise the repsonse to [0 1]
        GaborResp = SimplyNormalise((GaborResp).^2);
        % compute the scaling weight. Give higher priority to maps with low
        % median
        GaborResp = GaborResp./median(median(GaborResp));
        % take fft
        GaborResp = fftshift(fft2(GaborResp));
        % apply second order linear filtering (smoothing)
        GaborResp = B.*GaborResp;
        GaborResp = real(ifft2(ifftshift(GaborResp)));
        % add the map for a given level to the orientation channel map
        OriMap(:,:,ori) = OriMap(:,:,ori) + GaborResp;
    end
end
% imwrite(simplyNormalise(GaborFilters), 'GaborFilters.jpg');
end

function   SalMap = CombineOriMaps(OriMap, params, N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply weights to each orientation channel and add them to the overall
% saliency map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SalMap = zeros(N);
for ori = 1:params.NmbrOri
    m  = median(median(OriMap(:,:, ori)));
    % only count map if we have a non-zero median. Otherwise will be
    % dividing by 0 which is never fun.
    if m > 0
        ScalingWeight =max(max(OriMap(:,:, ori)))./m;
        OriMap(:,:, ori) = OriMap(:,:, ori).*ScalingWeight;
        SalMap =  SalMap + OriMap(:,:, ori)./params.NmbrOri;
    end
end
end



