function A = DownsamplingDisplacement(A, c)
% This function downsampling tensor A (m*n*p*3) by c and cast type as
% single for decreasing memory allcation, where p for phase

if isempty(A)
    return;
end

q = fix(c/2); % use the center point in block for resampling 
A = single(A);
A = pagetranspose(downsample(A,c,q));
A = pagetranspose(downsample(A,c,q));
end
