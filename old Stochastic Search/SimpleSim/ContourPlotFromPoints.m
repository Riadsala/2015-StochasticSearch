function ContourPlotFromPoints(data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transform data to right coord system
% [theta,r] = cart2pol(data(:,1), data(:,2));
% theta 
% clear data
% data = [theta, r];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% put into cart and smooth
% data = pol2cart(data(:,1), data(:,2));
min_edges = floor(min(data));
max_edges = ceil(max(data));
map = zeros(max_edges-min_edges+1);
data = ceil(data - repmat(min_edges, [size(data,1),1]))+1;
for t=1:size(data,1)
    map(data(t,1),data(t,2)) = map(data(t,1),data(t,2))+1;
end
 map = imfilter(map, fspecial('gaussian',[13 13], 1));
% figure
contour(map,20)
hold all
axis equal
plot([0, max_edges(2)-min_edges(2)]+1,[-min_edges(1),-min_edges(1)]+1, 'k-')
plot([-min_edges(2),-min_edges(2)]+1, [0, max_edges(1)-min_edges(1)]+1, 'k-')