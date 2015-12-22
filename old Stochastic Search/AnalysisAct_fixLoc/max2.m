function position = max2(Input)
% gives 2-d location of maximum value in input matrix
[m1 I1] = max(Input);
[m2 I2] = max(max(Input));
position = [I1(I2), I2];
end
