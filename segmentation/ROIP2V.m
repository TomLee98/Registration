function vroi = ROIP2V(rois, vsize)
%ROIP2V This function transform ROI representation to Ball representation
%by using linear indices in volume
%
% vroi = ROIP2V(rois, vsize)
%
% Input:
%   - rois:
%   - vsize:
% Output:
%   - vroi:

arguments
    rois (:,1) cell
    vsize (1,3) double {mustBeInteger, mustBePositive}
end

% transform xyz position to linear indices
vroi = cell(size(rois));
for n = 1:numel(rois)
    if ~isempty(rois{n})
        xyz_pts = [];
        for s = 1:size(rois{n}, 1)
            x0 = rois{n}(s, 1);
            y0 = rois{n}(s, 2);
            r = rois{n}(s, 3);
            z0 = rois{n}(s, 4);
            xy_pts = [];
            % calculate 1/4 circle
            for x = 0:r
                for y = 0:r
                    if x^2 + y^2 < r^2
                        xy_pts = [xy_pts; x,y]; %#ok<*AGROW>
                    end
                end
            end
            % using circle symmetry
            xy_pts = [xy_pts; [-xy_pts(:,1),xy_pts(:,2)]]; % x->-x
            xy_pts = [xy_pts; [xy_pts(:,1),-xy_pts(:,2)]]; % y->-y
            % translate the center
            xy_pts = round(xy_pts + [x0, y0]);
            % combine x,y,z
            xyz_pts = [xyz_pts; [xy_pts, z0*ones(size(xy_pts, 1), 1)]];
        end
        % transform to linear indices
        vroi{n} = sub2ind(vsize, xyz_pts(:,2), xyz_pts(:,1), xyz_pts(:,3));
    end
end

end

