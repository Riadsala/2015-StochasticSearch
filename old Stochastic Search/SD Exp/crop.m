function crop = crop(A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to crop a image matrix - remove excesse 0's on borders.
% Alasdair Clarke 20/6/06
% 
% for example
%     
%     0 0 0 0 0
%     0 1 0 0 0
%     0 1 1 1 0 
%     0 0 0 0 0
%     0 0 0 0 0 
%     
% would go to 
% 
%      1 0 0 
%      1 1 1
%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = size(A,2);
y = size(A,1);

%crop top
while A(1,:) == zeros(1, x)
    B=A(2:y,:);
    A=B;
    y = size(A,1);
end

%crop bottom
while A(y,:) == zeros(1, x)
    B=A(1:y-1,:);
    A=B;
    y = size(A,1);
end

%crop left
while A(:,1) == zeros(y, 1)
    B=A(:,2:x);
    A=B;
    x = size(A,2);
end

%crop right
while A(:,x) == zeros(y, 1)
    B=A(:,1:x-1);
    A=B;
    x = size(A,2);
end


crop=A;
end