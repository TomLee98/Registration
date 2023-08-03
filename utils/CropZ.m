function [C, L] = CropZ(A, ROI)
% CROPZ This function crop image stack A by using ROI on z direction
% Input:
%   - A: the gray image stack with z at the 4th dimension
%   - ROI: 1*n positive integer row array, with the selected z-plane indices
% Output:
%   - C: z-plane croped gray image stack
%   - L: number of selected z-planes at each volume

% Version 1.0.0
% Copyright (c) 2022-2023, Weihan Li

arguments
    A;
    ROI (1,:) double
end
if numel(size(A)) < 2
    error("Invalid Matrix A");
end
if ~issorted(ROI)
    error("Invalid time ROI, check the order");
end
if ~isempty(ROI)
    eval_str = "A(:,:,:,ROI)";
    add_str = repmat(",:",1,ndims(A)-4);
    % insert time cut at end
    eval_str = insertAfter(eval_str,eval_str.strlength-1,add_str.join(""));
    C = eval(eval_str);
else
    C = A;
end
L = numel(ROI);
end

