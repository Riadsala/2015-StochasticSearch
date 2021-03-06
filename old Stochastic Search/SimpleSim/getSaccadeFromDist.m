function [x y] = getSaccadeFromDist(l,t)

% read in data
shift = csvread('shift.csv');

t = ceil(t/5);
t = min(4, t);


switch l(1)
    case 1
        map = csvread(['Corner_t=' int2str(t) '.csv']);
    case 2
        map = csvread(['EdgeT_t=' int2str(t) '.csv']);
    case 3
        map = csvread(['EdgeL_t=' int2str(t) '.csv']);
    case 4
        map = csvread(['midSaccades_t=' int2str(t) '.csv']);
end

% scale to give a probability distribution
map = map./sum(map(:));
distribution = cumsum(map(:));
r=rand;
r = find(distribution>r,1);
[x y] = ind2sub(size(map), r);
x = (x + shift(1, (l(1)-1)*4+t))/5;
y = (y + shift(2, (l(1)-1)*4+t))/5;
x = VisualAngleToPixels(x);
y = VisualAngleToPixels(y);


switch l(1) %consider if location is corner
    case 1
        switch l(2)
            case 1
                x = -x;
                y =  y ;
            case 2
                x= -x ;
                y= -y;
            case 3;
                x=  x ;
                y=  y;
            case 4
                x=  x ;
                y= -y ;
                
        end
    case 2 % cosider horixontal edge
        switch l(2)
            case 1
                x = -x;
                y =  y ;
            case 2
                x =  x ;
                y =  y;
        end
    case 3 % cosider vertical edge
        switch l(2)
            case 1
                x =  x; 
                y =  y;
            case 2
               x =  x;
               y = -y;               
        end
    case 4
        x =  x;
        y =  y;
end


