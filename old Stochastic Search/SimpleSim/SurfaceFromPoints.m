function [map] =SurfaceFromPoints(data)
% close all

data(:,2) = ceil(mod(180/pi*data(:,2),360));
data(data(:,2)==0, 2) = 360;
data(:,1) = 5*data(:,1);
data = ceil(data);
max(data)
min(data)
map = zeros(100,360);

for t=1:size(data,1)
    map(data(t,1),data(t,2)) = map(data(t,1),data(t,2))+1;
end
map = imfilter(map, fspecial('gaussian',[51 51], 4), 'circular');
map(80:100,:) = 0;
% figure
% contour(map)
% hold all
% % axis([-5 20 -5 20])
% plot([0, max_edges(2)-min_edges(2)],[-min_edges(1),-min_edges(1)], 'k-')
% plot([-min_edges(2),-min_edges(2)], [0, max_edges(1)-min_edges(1)], 'k-')