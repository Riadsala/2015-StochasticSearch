function test
close all
t=1;
time=4;
dir = []
for temp = 1:250
    fix = [50,500];
    [l a b] = getBin(fix(t,1),fix(t,2));
    [x y] = getSaccadeFromDist3(l, time) ;
    x=x+fix(1);
    y=y+fix(2);
    fix = [fix; x y];    
    plot(fix(:,1), fix(:,2), 'g-');
    dir = [dir, cart2pol(x,y)];
    hold on
    axis([-50 1024, 0 1024]);
end

for temp = 1:250
    fix = [1000,500];
    [l a b] = getBin(fix(t,1),fix(t,2));  

    [x y] = getSaccadeFromDist3(l,time);  
    x=x+fix(1);
    y=y+fix(2);
    fix = [fix; x y];
    plot(fix(:,1), fix(:,2), 'rx-');
    %  axis([0 1024, 0 1024])
    hold on
end
for temp = 1:250
    fix = [500,50];
    [l a b] = getBin(fix(t,1),fix(t,2));    
    [x y] = getSaccadeFromDist3(l,time) ;
     x=x+fix(1);
    y=y+fix(2);
    fix = [fix; x y];
    plot(fix(:,1), fix(:,2), 'g-');
%     axis([0 1024, 0 1024])
    hold on
end
for temp = 1:250
    fix = [500,1002];
    [l a b] = getBin(fix(t,1),fix(t,2));    
    [x y] = getSaccadeFromDist3(l,time) ;     
    x=x+fix(1);
    y=y+fix(2);
    fix = [fix; x y];
    plot(fix(:,1), fix(:,2), 'b-');    
    hold on
end
end

