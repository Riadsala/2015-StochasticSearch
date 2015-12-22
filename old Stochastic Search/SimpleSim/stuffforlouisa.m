% stuff for Louisa
close all
clear all

% make a pretend foreground image
colourednoise = rand(500, 500, 3);
colourednoise(:,:,3) = 0;
foreground  = ones(750,750,3);
foreground(126:625, 126:625,:) = colourednoise;
foreground(325:425,325:425,:) = 1;
clear colourednoise
figure;
subplot(2,1,1)
imshow(foreground)

% make pretend background image
background = rand(750,750,3);
background(:,:,2:3) = 0;
subplot(2,1,2)
imshow(background)

% louisa can work out what the code below does.... :)

idx = find((foreground(:,:,1)<1).*(foreground(:,:,2)<1).*(foreground(:,:,3)<1));
[x,y] = ind2sub(size(foreground),idx);
compositeimage = background;
% shouldn't need a for loop for this, but for some odd reason I was getting
% an OUT OF MEMORY error on my computer. I'm pretty sure:
%compositeimage(x,y,:) = foreground(x,y,:);
% should work though. Odd.
for t=1:size(x,1)
    compositeimage(x(t),y(t),:) = foreground(x(t),y(t),:);
end


figure;
imshow(compositeimage)




