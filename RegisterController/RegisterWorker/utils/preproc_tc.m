function [vs, rois] = preproc_tc(vs, mfsize, dfsize, gfsize, ga)
%PREPROC_TC This funciton do image preprocess for better coarse registration
% Input:
%   - vs: m-by-n-by-p array, the image volume
%   - mfsize: 1-by-3 positive odd integer array, median filter size
%   - dfsize: 1-by-3 array, the dilate filter size, with [radius(pix), 
%             intensity threshold, volume size(pix) threshold]
%   - gfsize: 1-by-3 positive odd integer array, gaussian filter size
%   - ga: 1-by-1 nonnegtive double, 0 - 2, 1 as default, gamma
%         transformation coefficient(typical power transformation)
% Output:
%   - vs: m-by-n-by-p array, the processed image volume
%   - rois:n-by-6 double array, [x,y,z,w,h,d] as cuboid location in vs

arguments
    vs 
    mfsize  (1,3)   double {mustBeInteger, mustBeNonnegative} = [3,3,3]    % [x,y,z]
    dfsize  (1,3)   double {mustBeInteger} = [3,100,1000]               % [r,iTh,vTh]
    gfsize  (1,3)   double {mustBeInteger, mustBeNonnegative} = [3,3,1]    % [x,y,z]
    ga      (1,1)   double {mustBeInRange(ga, 0, 2)} = 1
end

% ============== remove possible salt and pepper noise =================
if all(mfsize > 0)
    vs = medfilt3(vs, mfsize, "replicate");
end

% ===== feature enhance: contrast balance by gamma transformation ======
vs = uint16(single(vs).^ga);

% =============== dilate filter for target enhancing ===================
if dfsize(1) > 0
    % create a mask with hard threshold
    if dfsize(2) >= 0
        vs_mask = (vs >= dfsize(2));
    else
        vs_u8 = Remap(vs, "uint8");
        T = graythresh(vs_u8);
        vs_mask = imbinarize(vs_u8, T); % otsu imbinarize
    end
    % dilate the image s.t. finely chopped connected
    vs_mask = imdilate(vs_mask, strel("sphere", dfsize(1)));

    if dfsize(3) >= 0
        % remove connected domains whose volume are lower than vTh
        vs_mask = bwareaopen(vs_mask, dfsize(3), 18);
    else
        % only keep the maximum volume object
        CC = bwconncomp(vs_mask, 18);
        S = regionprops3(CC, "Area");
        L = labelmatrix(CC);
        [~, maxloc] = max(S.Area);
        vs_mask = ismember(L, maxloc);
    end
    % modify the background to 0
    vs(~vs_mask) = 0;
end

if nargout == 2
    % require roi and dilate filter enabled
    if dfsize(1) > 0
        CC = bwconncomp(vs_mask, 18);
        rois = regionprops3(CC, "BoundingBox");
    else
        % whole volume
        rois = [1,1,1,size(vs, [2,1,3])];
    end
else
    rois = [];
end

% gaussian low pass filter for volume smooth, more robust
if all(gfsize > 0)
    vs = imgaussfilt3(vs, "FilterSize", gfsize, "Padding", "replicate");
end
end

