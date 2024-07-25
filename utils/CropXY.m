function [C, status] = CropXY(A, ROI)
% CROPXY This function crop image stack A by using ROI on xy plane
% Input:
%   - A: the gray image stack with X and Y at the 1st and 2nd dimension
%   - ROI: 2*2 positive integer matrix, with two points positions for rectangle
%     representation, which is [x1,y1; x2,y2]
% Output:
%   - C: xy croped gray image stack
%   - status: a struct with 'width' and 'height' of ROI (number of pixels)

% Version 1.0.0
% Copyright (c) 2022-2023, Weihan Li

arguments
    A
    ROI (2,2) double {mustBeInteger,mustBePositive} % (x1,y1;x2,y2)
end
if numel(size(A)) < 2
    error("Invalid Matrix A");
end
if ROI(2,1)>ROI(1,1) && ROI(2,2)>ROI(1,2) && ...
        (ROI(1,1)>=1 && ROI(2,1)<=size(A,2)) &&...
        (ROI(1,2)>=1 && ROI(2,2)<=size(A,1))
    eval_str = "A(ROI(1,2):ROI(2,2),ROI(1,1):ROI(2,1))";
    add_str = repmat(",:",1,ndims(A)-2);
    % insert slice cut at end
    eval_str = insertAfter(eval_str,eval_str.strlength-1,add_str.join(""));
    C = eval(eval_str);
else
    error("Invalid ROI");
end
status.width = ROI(2,1)-ROI(1,1)+1;
status.height = ROI(2,2)-ROI(1,2)+1;
end

