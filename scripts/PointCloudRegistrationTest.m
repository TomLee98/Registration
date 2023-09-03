%% This script for test point cloud registration method
% First block: load data and extract the nuclear position
% we use PNs as test sample
[opts,rt,movie,~] = loadfile("E:\si lab\Matlab Projects\CaImAn_Interface\data\test1_30.tif");

sc = squeeze(movie(:,:,1,:,:));
fc = squeeze(movie(:,:,2,:,:));

%%
% col-1: affine3D object, col-2: 4D displacement field tensor
tforms_affine = cell(opts.frames, 1);
tforms_affine{1} = affinetform3d(eye(4));
tforms_nonrigid = cell(opts.frames, 1);
tforms_nonrigid{1} = zeros([size(sc,[1,2,3]), 3]);
rmses = zeros(opts.frames, 1);
%% generate the key frame registration
Rs = [opts.xRes, opts.yRes, opts.zRes];
RA = imref3d(size(sc(:,:,:,1)), Rs(1), Rs(2), Rs(3));
pts_init = vol2pts(sc(:,:,:,1), Rs);
fix_init = sc(:,:,:,1);
pts_init_ds = pcdownsample(pts_init, "gridAverage",2);
fix_init_ds = DownSampling(sc(:,:,:,1), 0.25);
RA_ds = imref3d(size(fix_init_ds), Rs(1)/0.25, Rs(2)/0.25, Rs(3));
[optimizer, metric] = imregconfig("monomodal");
% optimizer.MaximumIterations = 100;
optimizer.MinimumStepLength = 1e-7;
optimizer.MaximumStepLength = 0.0625;
optimizer.RelaxationFactor = 0.5;

for k = 1:opts.frames-1
    % generate points cloud
    fixvol = sc(:,:,:,k);
    movol = sc(:,:,:,k+1);
    movol_ds = DownSampling(movol, 0.25);
    MOV_BKG = prctile(movol, 10, "all");

    pts_fix = vol2pts(fixvol, Rs);
    pts_mov = vol2pts(movol, Rs);

    % using cloud registration cpd -> affine registration
    % downsampling, rough registration, make sure the structured space
    pts_fix_ds = pcdownsample(pts_fix, "gridAverage",2);    % 2 um for radius, the sampling theorem
    pts_mov_ds = pcdownsample(pts_mov, "gridAverage",2);

    % first estimation for key k -> key k+1 points cloud registration
    [tform_affine, ~, rmses(k+1)] = pcregistercpd(pts_mov_ds, pts_fix_ds, ...
        "Transform","Rigid","OutlierRatio",0.1,"MaxIterations",100,"Tolerance",1e-7);
    tforms_affine{k+1} = affinetform3d();

    % pcshowpair(pts_fix_ds, pctransform(pts_mov_ds, tform_affine),"MarkerSize",30);
    
    % update the fixed -> k+1 transformation estimation
    tforms_affine{k+1}.A = tform_affine.A*tforms_affine{k}.A;

    % second update tforms again to avoid accumulation of errors
    tforms_affine{k+1} = imregtform(movol_ds, RA_ds, fix_init_ds, RA_ds, "affine", ...
        optimizer, metric, "InitialTransformation", tforms_affine{k+1});

    % pcshowpair(pts_init_ds, pctransform(pts_mov_ds, tforms_affine{k+1}),"MarkerSize",30);

    % update the fixed -> k+1 fine tuned nonrigid estimation
    reg_affined = imwarp(movol, RA, tforms_affine{k+1},"linear",...
        "OutputView",RA,"FillValues",MOV_BKG);

    [tforms_nonrigid{k+1}, ~] = imregdemons(reg_affined, fix_init, [100,50,25],...
            "AccumulatedFieldSmoothing",1,"DisplayWaitbar",false);

    fprintf("progress: %.1f %%\n", k/(opts.frames-1)*100);
end

%% view total results
figure;
for n = 1:opts.frames
    reg_affined = imwarp(sc(:,:,:,n),RA,tforms_affine{n},"linear","OutputView",RA,...
        "FillValues",MOV_BKG);
    reg_nonrigid = imwarp(reg_affined, tforms_nonrigid{n},"cubic",...
        "FillValues",MOV_BKG,"SmoothEdges",true);
    for k = 1:size(reg_nonrigid, 3)
        imshowpair(sc(:,:,k,1),reg_affined(:,:,k),"falsecolor");
        title(sprintf("volume: %d, slice: %d", n, k),"FontSize",15);
        pause(0.05);
    end
end


%%
function J = remap(I, type, LUT)
% This function remap the image I intensity to new type by defined LUT and
% it will fully map the target gray scale
arguments
    I {mustBeNonnegative};
    type (1,1) string {mustBeMember(type, ["uint8","uint16","uint32","single","double"])} = "uint8";
    LUT (1,2) double {mustBeNonnegative} = [0, inf];
end

if all(LUT == [0, inf])
    % auto scale the data by min and max
    LUT = [min(I,[],"all"), max(I,[],"all")];
else
    if LUT(1) < min(I,[],"all") || LUT(2) > max(I,[],"all")
        warning("Auto shrink the LUT to fit the gray scale.");
        LUT(1) = max(LUT(1), min(I,[],"all"));
        LUT(2) = min(LUT(2), max(I,[],"all"));
    end
end

validateattributes(LUT,{'numeric'},{'increasing', 'row'});

% remap the intensity
J = double(I - LUT(1))/double(LUT(2)- LUT(1));

if ismember(type, ["uint8","uint16","uint32"])
    J = cast(J*double(intmax(type)), type);
elseif type == "single"
    J = single(J);
end

end

function pts = vol2pts(vol, res, minobj, rlut)
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
vol_u8 = remap(vol, "uint8");
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