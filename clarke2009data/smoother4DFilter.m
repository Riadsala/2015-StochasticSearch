function F = smoother4DFilter(N, sigma)

for w = 1:N
    for x = 1:N
        for y = 1:N
            for z = 1:N              
                F(w,x,y,z) = mvnpdf([w x y z], (N+1)/2 * ones(1,4), sigma*eye(4));


            end
        end
    end
end