function tform = imregmc(moving, fixed, Rf)
% IMREGMMT1: This function use 1 order moment (intensity mean) for 
% translation shift estimation
% Input:
%   - moving:
%   - fixed:
%   - Rf:
% Output:
%   - tform:
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

mc_ref_x = sum((1:nx).*sum(fixed, 1))/mfixed;
mc_ref_y = sum((1:ny)'.*sum(fixed, 2))/mfixed;

mc_moving_x = sum((1:nx).*sum(moving, 1))/mmoving;
mc_moving_y = sum((1:ny)'.*sum(moving, 2))/mmoving;

% translate the pixels to image coordinate system
tf_x = (mc_moving_x - mc_ref_x)*Rf.PixelExtentInWorldX;
tf_y = (mc_moving_y - mc_ref_y)*Rf.PixelExtentInWorldY;

if isMATLABReleaseOlderThan("R2022b")
    T = [[eye(2);[tf_x, tf_y]],[0;0;1]];
    tform = affine2d(T);
else
    tform = transltform2d(tf_x, tf_y);
end
end