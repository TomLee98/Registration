function [C, L] = CropT(A, ROI)
% CROPT This function crop image stack A by using ROI on time sequence
% Input:
%   - A: the gray image stack with time at the last dimension
%   - ROI: 1*n positive integer row array, with the selected frame indices
% Output:
%   - C: time croped gray image stack
%   - L: number of selected frames

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
    eval_str = "A(ROI)";
    add_str = repmat(":,",1,ndims(A)-1);
    % insert time cut at end
    eval_str = insertAfter(eval_str,2,add_str.join(""));
    C = eval(eval_str);
else
    C = A;
end
L = numel(ROI);
end

