function vs = preproc_oc(vs, mfsize, gfsize)
%PREPROC_OC This funciton do image preprocess for better coarse registration
% Input:
%   - vs: m-by-n-by-p array, the image volume
%   - mfsize: 1-by-3 positive odd integer array, median filter size
%   - gfsize: 1-by-3 positive odd integer array, gaussian filter size
% Output:
%   - vs: m-by-n-by-p array, the processed image volume

arguments
    vs
    mfsize  (1,3)   double {mustBeInteger, mustBeNonnegative} = [3,3,3]     % [x,y,z]
    gfsize  (1,3)   double {mustBeInteger, mustBeNonnegative} = [3,3,1]     % [x,y,z]
end

% ============== remove possible salt and pepper noise =================
if all(mfsize > 0)
    vs = medfilt3(vs, mfsize, "replicate");
end

% gaussian low pass filter for volume smooth, more robust
if all(gfsize > 0)
    vs = imgaussfilt3(vs, "FilterSize", gfsize, "Padding", "replicate");
end

end