function [vs, rois] = preproc_tc(vs, dfsize, mfsize, gfsize, gamma, roiref)
%PREPROC_TC This funciton do image preprocess for better coarse registration
% Input:
%   - vs: m-by-n-by-p array, the image volume
%   - mfsize: 1-by-3 positive odd integer array, median filter size
%   - dfsize: 1-by-3 array, the dilate filter size, with [radius(pix), 
%             intensity threshold, volume size(pix) threshold]
%   - gfsize: 1-by-3 positive odd integer array, gaussian filter size
%   - gamma: 1-by-1 nonnegtive double, 0 - 2, 1 as default, gamma
%         transformation coefficient(typical power transformation)
%   - roiref: k-by-3 array, with [w,h,d] for objects in reference scene
% Output:
%   - vs: m-by-n-by-p array, the processed image volume
%   - rois:n-by-6 double array, [x,y,z,w,h,d] as cuboid location in vs

arguments
    vs 
    dfsize  (1,3)   double {mustBeInteger} = [3,100,1000]                   % [r,iTh,vTh]
    mfsize  (1,3)   double {mustBeInteger, mustBeNonnegative} = [3,3,3]     % [x,y,z]
    gfsize  (1,3)   double {mustBeInteger, mustBeNonnegative} = [3,3,1]     % [x,y,z]
    gamma   (1,1)   double {mustBeInRange(gamma, 0, 2)} = 1
    roiref  (:,3)   double {mustBePositive, mustBeInteger} = double.empty(0,3)
end

% ============== remove possible salt and pepper noise =================
if all(mfsize > 0)
    vs = medfilt3(vs, mfsize, "replicate");
end

% ===== feature enhance: contrast balance by gamma transformation ======
vs = uint16(single(vs).^gamma);

% =============== dilate filter for target enhancing ===================
if ~isempty(roiref)
    % calculate the best dilatefilter for fitting the reference
    % optimize the vs_mask until the distance goes minimum
    % optimizer: genetic algorithm
    drm = max(round((3/(4*pi)*prod(roiref(end, :)))^(1/3)), 5);
    lb = [1, 1/32*pi*prod(roiref(end,:))];
    ub = [drm,  1/4*pi*prod(roiref(end,:))];
    vs_min = min(vs,[],"all"); vs_max = max(vs, [], "all");
    vs_u8 = Remap(vs, "uint8");
    T = graythresh(vs_u8);
    vm = imbinarize(vs_u8, T);     % otsu imbinarize
    opf = @(x)whdist(x, roiref, vm);
    opts = optimoptions("ga", "PopulationSize",20, "MaxGenerations",5);
    dfsize = ga(opf, 2, [],[],[],[], lb,ub, [], [1,2], opts);
    dfsize = double([dfsize(1), T*(vs_max-vs_min)+vs_min, dfsize(2)]);
end

if dfsize(1) > 0
    % create a mask with hard threshold
    if dfsize(2) >= 0
        vs_mask = (vs >= dfsize(2));
    else
        vs_u8 = Remap(vs, "uint8");
        T = graythresh(vs_u8);
        vs_mask = imbinarize(vs_u8, T);     % otsu imbinarize
    end
    % dilate the image s.t. finely chopped connected
    vs_mask = imdilate(vs_mask, strel("sphere", dfsize(1)));

    if dfsize(3) >= 0
        % remove connected domains whose volume are lower than vTh
        vs_mask = bwareaopen(vs_mask, dfsize(3), 18);
    else
        % only keep the maximum volume object
        CC = bwconncomp(vs_mask, 18);
        S = regionprops3(CC, "Volume");
        L = labelmatrix(CC);
        [~, maxloc] = max(S.Volume);
        vs_mask = ismember(L, maxloc);
    end
    % modify the background to 0
    vs(~vs_mask) = 0;
end

if nargout == 2
    % require roi and dilate filter enabled
    if dfsize(1) > 0
        CC = bwconncomp(vs_mask, 18);
        rois = regionprops3(CC, "BoundingBox", "Volume");
        if ~isempty(rois)
            % sort the rois by volume
            [~, vmidx] = sort(rois.Volume, 1, "descend");
            rois = rois.BoundingBox(vmidx, :);
        else
            rois = double.empty(0, 6);
        end
    else
        % whole volume
        rois = [1,1,1,size(vs, [2,1,3])];
    end
else
    rois = double.empty(0, 6);
end

% gaussian low pass filter for volume smooth, more robust
if all(gfsize > 0)
    vs = imgaussfilt3(vs, "FilterSize", gfsize, "Padding", "replicate");
end
end

function f = whdist(x, u0, v)
% This function apply x as dilatefilter on v, which could return a border
% matrix compare with u0
% note that x(1) must be positive

% binarize the volume by global hard threshold
vm = (v >= x(2));

% dilate binarized volume
vm = imdilate(vm, strel("sphere", x(1)));

% remove connected domains whose volume are too small
vm = bwareaopen(vm, x(2), 18);

% calculate the object borders 
CC = bwconncomp(vm, 18);
rois = regionprops3(CC, "BoundingBox", "Volume");
if ~isempty(rois)
    % sort the rois by volume
    [~, vmidx] = sort(rois.Volume, 1, "descend");
    rois = rois(vmidx, :);
    u = rois.BoundingBox(:, [5,4,6]);
    nobj = size(u, 1); nobj0 = size(u0, 1);
    if nobj < nobj0
        % add zeros behind u
        u = [u; zeros(nobj0-nobj, 3)]; 
    elseif nobj > nobj0
        u0 = [u0; zeros(nobj-nobj0, 3)]; 
    end
else
    u = zeros(size(u0));
end

f = sum((u-u0).^2, "all");

end