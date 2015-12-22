function params = getParams
maxRes = 8;
params.start_l = 1;
params.NumberOfRes = maxRes -params.start_l;
params.orientations = [0 22.5, 45, 67.5, 90, 112.5, 135, 157.5];

params.NmbrOri = size(params.orientations,2);
params.Gabor_stddev =1.9;
params.Gabor_r = 1.8;

params.EDPradius = -0.001;
params.nmbrOfFixations = 20;
params.IOR_mask_size = 30;
params.IOR_decaytime = 10;
params.FixationRadius = 1; %given in degrees of visual angles

params.FixationRadius = VisualAngleToPixels(params.FixationRadius);

end


function pixeldist = VisualAngleToPixels(VisualAngle)
pixelsize = 0.255; % in mm
distfromdisplay = 0.87; % in meters
ratio = tand(VisualAngle);
dist = ratio*distfromdisplay;
pixeldist = dist/(pixelsize/1000);
end