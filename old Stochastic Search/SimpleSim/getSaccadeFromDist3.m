function [x y] = getSaccadeFromDist3(l,t)

% read in data

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

r=ceil(size(map,1)*rand);
saccade = map(r,:);

[y x] = pol2cart(saccade(2), saccade(1));
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


