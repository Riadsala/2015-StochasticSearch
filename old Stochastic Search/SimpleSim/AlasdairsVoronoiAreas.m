function Area = AlasdairsVoronoiAreas(fix)


%cheat and just assign pixels to nearest fixation
% fix = 1024*rand(10,2);
n = size(fix,1);
X=repmat(1:1024,[1024 1]);
Y=X';
Dmap = zeros(1024,1024,n);
for f=1:n
    Dmap(:,:,f) = sqrt((X-fix(f,1)).^2+(Y-fix(f,2)).^2);
end


Assignments = zeros(1024);
Area = zeros(n,1);
for f=1:n
    D = min(Dmap(:,:,1:f), [], 3);    
    for t=1:f
         Assignments(D==Dmap(:,:,f)) = f;
        Area(f) = max(Area(f), sum(sum(D==Dmap(:,:,t))));
    end
%     figure
% imshow(Assignments,[])
end

% figure;
% plot(Area)
