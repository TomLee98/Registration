function tform = imregmc3(moving, fixed, Rf, cargs)
% IMREGMMT1: This function use 1 order moment (intensity mean) for 
% translation shift estimation
% Input:
%   - moving: m-by-n-by-p uint16 matrix, the moving image
%   - fixed: m-by-n-by-p uint16 matrix, the reference image
%   - Rf: 1-by-1 imref3d object, as reference coordinate system
% Output:
%   - tform: 1-by-1 transltform3d object

arguments
    moving  (:,:,:) uint16
    fixed   (:,:,:) uint16
    Rf      (1,1)   imref3d
    cargs   (1,1)   struct = struct("Filter", "median", "VT", 1000, "Radius", 3)
end

if any(size(moving)~=size(fixed))
    throw(MException("imregmmt1:imageSizeNotMatch", ...
        "Moving and fixed image must have the same size."));
end

[ny, nx, nz] = size(fixed);

% functional channel, with large local fluorescence mutation
% use binary mask for robust coarse estimation
fixed_omsk = isoutlier(single(fixed), cargs.Filter);
moving_omsk = isoutlier(single(moving), cargs.Filter);

% volume open operation to remove possible local little structure
fixed_omsk = bwareaopen(fixed_omsk, cargs.VT, 26);
moving_omsk = bwareaopen(moving_omsk, cargs.VT, 26);

% erode for better true foreground estimation
fixed_omsk = imerode(fixed_omsk, strel("sphere", cargs.Radius));
moving_omsk = imerode(moving_omsk, strel("sphere", cargs.Radius));

% calculate total volumes
mfixed = sum(fixed_omsk, "all");
mmoving = sum(moving_omsk, "all");

if abs(mmoving - mfixed)/mfixed > 0.5
    warning("imregmc3:badForegroundEstimation", "Mass center processing did not " + ...
        "obtain a similar size. The resulting registration may not be ideal.");
end

% calculate the mass center
mc_ref_x = sum((1:nx).*sum(fixed_omsk, [1,3]))/(mfixed+eps);
mc_ref_y = sum((1:ny)'.*sum(fixed_omsk, [2,3]))/(mfixed+eps);
mc_ref_z = sum(reshape(1:nz,1,1,[]).*sum(fixed_omsk,[1,2]))/(mfixed+eps);

mc_moving_x = sum((1:nx).*sum(moving_omsk, [1,3]))/(mmoving+eps);
mc_moving_y = sum((1:ny)'.*sum(moving_omsk, [2,3]))/(mmoving+eps);
mc_moving_z = sum(reshape(1:nz,1,1,[]).*sum(moving_omsk,[1,2]))/(mmoving+eps);

% translate the pixels to image coordinate system
tf_x = (mc_ref_x - mc_moving_x)*Rf.PixelExtentInWorldX;
tf_y = (mc_ref_y - mc_moving_y)*Rf.PixelExtentInWorldY;
tf_z = (mc_ref_z - mc_moving_z)*Rf.PixelExtentInWorldZ;

tform = transltform3d(tf_x, tf_y, tf_z);
end
