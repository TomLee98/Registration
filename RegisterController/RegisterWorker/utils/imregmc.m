function tform = imregmc(moving, fixed, Rf)
% IMREGMMT1: This function use 1 order moment (intensity mean) for 
% translation shift estimation
% Input:
%   - moving: m-by-n uint16 matrix, the moving image
%   - fixed: m-by-n uint16 matrix, the reference image
%   - Rf: 1-by-1 imref2d object, as reference coordinate system
% Output:
%   - tform: 1-by-1 affine2d(MATLAB<R2022b) or transltform2d object
%

arguments
    moving  (:,:)   uint16
    fixed   (:,:)   uint16
    Rf      (1,1)   imref2d
end

if any(size(moving)~=size(fixed))
    throw(MException("imregmmt1:imageSizeNotMatch", ...
        "Moving and fixed image must have the same size."));
end

[ny, nx] = size(fixed);
mfixed = sum(fixed, "all");
mmoving = sum(moving, "all");

mc_ref_x = sum((1:nx).*sum(fixed, 1))/(mfixed+eps);
mc_ref_y = sum((1:ny)'.*sum(fixed, 2))/(mfixed+eps);

mc_moving_x = sum((1:nx).*sum(moving, 1))/(mmoving+eps);
mc_moving_y = sum((1:ny)'.*sum(moving, 2))/(mmoving+eps);

% translate the pixels to image coordinate system
tf_x = (mc_ref_x - mc_moving_x)*Rf.PixelExtentInWorldX;
tf_y = (mc_ref_y - mc_moving_y)*Rf.PixelExtentInWorldY;

tform = transltform2d(tf_x, tf_y);
end