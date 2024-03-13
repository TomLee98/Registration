function mask = recoverMask(mask, sz)
%RECOVERMASK This function recover the sparse matrix to 3d array

arguments
    mask    (:,:)
    sz      (1,3)   double  {mustBePositive, mustBeInteger}
end

if ~issparse(mask)
    throw(MException("recoverMask:invalidMaskFormat", ...
        "Input mask must be sparse matrix."));
end
if size(mask, 1) ~= prod(sz)
   throw(MException("recoverMask:invalidMaskSize", ...
        "Invalid mask size."))
end

% reconstruct the mask volume
v_mask = zeros(sz, "like", mask);
for k = 1:size(mask, 2)
    vidx = full(mask(:, k));
    v_mask(vidx>0) = k;
end

mask = v_mask;

end

