function [l, a,b] = getBin(x,y)
a = 0;
b = 0;
if (x>0)&&(x<205)
    a=1;
elseif (x>=205)&&(x<410)
    a=2;
elseif (x>=410)&&(x<615)
    a=3;
elseif (x>=615)&&(x<820)
    a=4;
elseif (x>=820)&&(x<1025)
    a=5;
end
if (y>0)&&(y<205)
    b=1;
elseif (y>=205)&&(y<410)
    b=2;
elseif (y>=410)&&(y<615)
    b=3;
elseif (y>=615)&&(y<820)
    b=4;
elseif (y>=820)&&(y<1025)
    b=5;
end
if (a==1)&&(b==1)
    l=[1 1];
elseif (a==1)&&(b==5)
    l=[1 2];
elseif (a==5)&&(b==1)
    l=[1 3];
elseif (a==5)&&(b==5)
    l=[1 4];
elseif ((b==2)+(b==3)+(b==3))&&(a==1)
    l=[2 1];
elseif ((b==2)+(b==3)+(b==3))&&(a==5)
    l=[2 2];
elseif ((a==2)+(a==3)+(a==3))&&(b==1)
    l=[3 1];
elseif ((a==2)+(a==3)+(a==3))&&(b==5)
    l=[3 2];
else
    l = [4 1];
end
end