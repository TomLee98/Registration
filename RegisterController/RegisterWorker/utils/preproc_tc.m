function vs = preproc_tc(vs, mfsize, ofsize, gfsize, ga)
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
    mfsize  (1,3)   double {mustBeInteger, mustBePositive} = [3,3,3]
    ofsize  (1,2)   double {mustBeInteger, mustBePositive} = [3,2]
    gfsize  (1,3)   double {mustBeInteger, mustBePositive} = [3,3,1]
    ga      (1,1)   double {mustBeInRange(ga, 0, 2)} = 1
end

% remove possible salt and pepper noise
vs = medfilt3(vs, mfsize, "replicate");

% feature enhance: contrast balance by gamma transformation
vs = uint16(single(vs).^ga);

% open operation to remove strong bright noise
vs = imopen(vs, offsetstrel("ball", ofsize(1), ofsize(2)));

% comment this, too memory allocated
% imhistmatch for photobleaching recovery
% vs = imhistmatchn(vs_, refv_, hn_);   

% gaussian low pass filter for volume smooth, more robust
vs = imgaussfilt3(vs, "FilterSize", gfsize, "Padding", "replicate");
end

