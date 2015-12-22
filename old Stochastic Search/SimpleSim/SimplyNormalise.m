function output = SimplyNormalise(I)

output = I - min(I(:));
if max(output(:))~=0
    output = output./max(output(:));
end
end