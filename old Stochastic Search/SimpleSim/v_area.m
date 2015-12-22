 function output=v_area(X)


[v c] =voronoin(X);
for j=1:length(c)
A(j) = polyarea(v(c{j},1), v(c{j},2));
end
output = A(isfinite(A));

