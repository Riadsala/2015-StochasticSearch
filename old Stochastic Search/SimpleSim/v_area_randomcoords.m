subplot(3,3,9)
cells = []
 for t=1:20
 
X = 1024*rand(32,2);
[v c] =voronoin(X);
for j=1:length(c)
A(j) = polyarea(v(c{j},1), v(c{j},2));
end
output = A(isfinite(A));
cells = [cells, output];
 end
 median(cells)
 cells(cells>500000)=[];
  h=hist(cells,40);
%     h=hist(log(V_Areas), 0.25:0.5:19.75);    
     h = h./size(cells,2);
     bar(625:1250:50000, h, 1,'g')    