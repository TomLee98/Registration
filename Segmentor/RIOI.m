function index = RIOI(rois, bnds)
%RIOI This function allocate ROIs indices to match subscreen boundaries
% Input:
%   - rois:
%   - bnds:
% Output:
%   - index:

arguments
    rois (:,1) cell
    bnds (:,4) double {mustBeInteger, mustBePositive}
end

index = ones(numel(rois), 1);
A = cell(size(bnds,1), 1);
for s = 1:size(bnds, 1)
    A{s} = [bnds(s,1), bnds(s,1)+bnds(s,3)-1; ...
            bnds(s,2), bnds(s,2)+bnds(s,4)-1];
end


for n = 1:numel(rois)
    B = [min(rois{n}(:,1)-rois{n}(:,3)), max(rois{n}(:,1)+rois{n}(:,3)); ...
         min(rois{n}(:,2)-rois{n}(:,3)), max(rois{n}(:,2)+rois{n}(:,3))];
    for s = 1:size(bnds, 1)
        C = A{s} - B;
        if (C(1,1)<=0)&&(C(2,1)<=0)&&(C(1,2)>=0)&&(C(2,2)>=0)
            index(n) = s;
            break;
        end
    end
end

end

