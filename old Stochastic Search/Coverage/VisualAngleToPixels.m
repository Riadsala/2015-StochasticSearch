function pixeldist = VisualAngleToPixels(VisualAngle)
pixelsize = 0.255; % in mm
distfromdisplay = 0.87; % in meters
ratio = tand(VisualAngle);
dist = ratio*distfromdisplay;
pixeldist = dist/(pixelsize/1000);
end