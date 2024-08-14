function [R, W, H] = ROI2Rectangle(ROI)
%ROI2RECTANGLE This function return the coordinate in UIAxes, real number
% Input:
%   - ROI: images.roi, could be images.roi.Circle, images.roi.Freehand,
%     images.roi.Polygon, images.roi.Rectangle
% Output:
%   - R: the 2*2 ROI matrix, pixels
%   - W: the width of ROI, pixels
%   - H: the height of ROI, pixels

% Version 1.0.0
% Copyright (c) 2022-2023, Weihan Li

arguments
    ROI (1,1) {mustBeA(ROI,{'images.roi.Circle','images.roi.Freehand',...
        'images.roi.Polygon','images.roi.Rectangle'})};
end

switch class(ROI)
    case "images.roi.Circle"
        R = [ROI.Center-ROI.Radius;ROI.Center+ROI.Radius];
        W = 2*ROI.Radius;
        H = 2*ROI.Radius;
    case "images.roi.Rectangle"
        R = [ROI.Vertices(1,:);ROI.Vertices(3,:)];
        W = ROI.Vertices(3,1)-ROI.Vertices(1,1);
        H = ROI.Vertices(3,2)-ROI.Vertices(1,2);
    otherwise
        R = [min(ROI.Position(:,1)),min(ROI.Position(:,2));...
            max(ROI.Position(:,1)),max(ROI.Position(:,2))];
        W = max(ROI.Position(:,1))-min(ROI.Position(:,1));
        H = max(ROI.Position(:,2))-min(ROI.Position(:,2));
end
end

