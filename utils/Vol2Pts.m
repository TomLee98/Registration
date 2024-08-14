function pts = Vol2Pts(vol, res, minobj, rlut)
% This function does pre-processing and transformes voxels to point cloud
arguments
    vol (:,:,:);
    res (1,3) double {mustBePositive};
    minobj (1,1) double {mustBePositive, mustBeInteger} = 1000;
    rlut (1,2) double {mustBeInRange(rlut, 0,100)} = [99.99, 97];
end

% ========================== preprocessing ===========================
% limit the over exposure voxels
vol(vol>prctile(vol,rlut(1),"all")) = prctile(vol,rlut(1),"all");

% remove background baseline
vol = vol - median(vol,"all");

% using otsu algorithm for calculating nuclear pixels
vol_u8 = Remap(vol, "uint8");
clear vol;

% quantile
vol_u8(vol_u8<prctile(vol_u8,rlut(2),"all")) = 0;

% remove little object on the foreground
vol_bw = bwareaopen(imbinarize(vol_u8,"global"), minobj, 26);
clear vol_u8;

% transform bw to xyz
vol_pts_idx = find(vol_bw);
[vps_y, vps_x, vps_z] = ind2sub(size(vol_bw), vol_pts_idx);
clear vol_bw;

% note that coordinate resolution need to be introduced here
% so the point with world scale: um as unit
pts = [vps_x, vps_y, vps_z].*res;

% generate two point clouds
pts = pointCloud(pts);
end