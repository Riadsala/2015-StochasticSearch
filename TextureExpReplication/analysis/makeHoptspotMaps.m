fix = round(csvread('allFixations.txt',1));

ctr = 0;
for person = [2 4 5]
    ctr = ctr + 1;
    pfix = fix(fix(:,1)==person,:);
    map = zeros(1024, 1024);
    for f = 1:length(pfix)
        map(pfix(f,3), pfix(f,2)) = map(pfix(f), pfix(f))+1;
    end
    B = fspecial('gaussian', [101 101], 19);
    map = imfilter(map, B);
    subplot(1,3,ctr)
    imshow(map,[])
    colormap('hot')
    size(map)
end
set(gcf, 'Color', 'w');
export_fig('hotspotmaps.pdf')