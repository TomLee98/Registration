function muinfo = MutualInfo(R,F,method)
% MUTUALINFO This function calculate the mutual information between image R
% and image F, which required grayscale or RGB image
% Input:
%   - R: the image R, type must be numeric, which could be RGB image,
%           dimension free but equal to F
%   - F: the image R, type must be numeric, which could be RGB image,
%           dimension free but equal to R
%   - method: the mutual information calculation method, could be "NMI" or
%           "MI", "NMI" as default
% Output:
%   - muinfo: the mutual information between image R and image F

% Version 1.0.0
% Copyright (c) 2022-2023, Weihan Li

arguments
    R {mustBeNumeric}
    F {mustBeNumeric}
    method (1,1) string {mustBeMember(method, ["MI","NMI"])} = "NMI";
end

assert(ndims(R)==ndims(F),"MutualInfo:dimensionNotMatch","Dimension not match.");
if ndims(R) == 3
    % transform to 2D problem
    if size(F,3) == 3
        F = rgb2gray(F);
    end
    if size(R,3) == 3
        R = rgb2gray(R);
    end
else
    %
end

% calculate the mutual information

% rescale the image to 10-bits image by min/max linear transformation
R_min = min(R,[],"all"); R_max = max(R,[],"all");
F_min = min(F,[],"all"); F_max = max(F,[],"all");
R = uint16(double(R - R_min)/(double(R_max - R_min) + eps)*1023 + 1);
F = uint16(double(F - F_min)/(double(F_max - F_min) + eps)*1023 + 1);

% calculate the empirical joint probability distribution
P = histcounts2(R, F, "BinWidth",1,"Normalization","probability");

Marg_A = sum(P); %sum at each column for edge probability distribution
Marg_B = sum(P,2); %sum at each row for edge probability distribution

% calculate the information entropy
H_A = sum(-Marg_A.*log2(Marg_A + eps));
H_B = sum(-Marg_B.*log2(Marg_B + eps));

% calculate the joint information entropy
H_AB = sum(-P.*log2(P + eps),"all");

switch method
    case "MI"
        muinfo = H_A + H_B - H_AB;
    case "NMI"
        muinfo = 0.5*(H_A + H_B)/H_AB;
    otherwise
end

end