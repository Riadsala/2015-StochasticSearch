function output = PixelsToVisualAngle(input)

% input is a distance given in pixels

pixelsize = 0.255; % in mm
distfromdisplay = 0.87; % in meters

saccadelength = (pixelsize/1000)*input;

output = atand(saccadelength./distfromdisplay);

end