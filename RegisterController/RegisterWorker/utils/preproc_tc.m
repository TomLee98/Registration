function vs = preproc_tc(vs, mfsize, dfsize, gfsize, ga)
%PREPROC_TC This funciton do image preprocess for better coarse registration
% Input:
%   - vs:
%   - mfsize:
%   - ofsize:
%   - gfsize:
%   - ga:
% Output:
%   - vs:
arguments
    vs 
    mfsize  (1,3)   double {mustBeInteger, mustBePositive} = [3,3,3]    % [x,y,z]
    dfsize  (1,3)   double {mustBeInteger} = [3,100,1000]               % [r,iTh,vTh]
    gfsize  (1,3)   double {mustBeInteger, mustBePositive} = [3,3,1]    % [x,y,z]
    ga      (1,1)   double {mustBeInRange(ga, 0, 2)} = 1
end

% ============== remove possible salt and pepper noise =================
vs = medfilt3(vs, mfsize, "replicate");

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

    end
    % modify the background to 0
    vs(~vs_mask) = 0;
end

% gaussian low pass filter for volume smooth, more robust
vs = imgaussfilt3(vs, "FilterSize", gfsize, "Padding", "replicate");
end

