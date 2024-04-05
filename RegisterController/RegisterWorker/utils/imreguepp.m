function tf_est = imreguepp(moving, fixed, rsFixed, tfType)
%IMREGUEPP This function use point cloud registration for motion correction
% Input:
%   - moving:
%   - fixed:
%   - rsFixed:
%   - tfType:
% Output:
%   - tf_est:
%
%   see also: pointCloud, pcregistercpd, bwulterode, imbinarize

% Create binarized representation
moving_u8 = Remap(moving, "uint8");
bw_moving = imbinarize(moving_u8, "global");
bw_moving_ue = bwulterode(bw_moving, "euclidean", 26);
moving_ity = reshape(moving(bw_moving_ue), [], 1);

fixed_u8 = Remap(fixed, "uint8");
bw_fixed = imbinarize(fixed_u8, "global");
bw_fixed_ue = bwulterode(bw_fixed, "euclidean", 26);
fixed_ity = reshape(fixed(bw_fixed_ue), [], 1);

% Generate point cloud and matching
xyz_moving = find(bw_moving_ue);
[y, x, z] = ind2sub(size(bw_moving), xyz_moving);
pt_moving = pointCloud([x,y,z].*rsFixed, "Intensity",moving_ity);

xyz_fixed = find(bw_fixed_ue);
[y, x, z] = ind2sub(size(bw_fixed), xyz_fixed);
pt_fixed = pointCloud([x,y,z].*rsFixed, "Intensity",fixed_ity);

switch tfType
    case {'translation', 'rigid'}
        tf_est = pcregistercpd(pt_moving, pt_fixed, Transform="Rigid", ...
            OutlierRatio=0.2, MaxIterations=15);    % rigidtform3d
        if "translation" == tfType
            % transformation cast
            tf_est = transltform3d(tf_est.Translation); % transltform3d
        end
    case "affine"
        tf_est = pcregistercpd(pt_moving, pt_fixed, Transform="Affine", ...
            MaxIterations=30);
end

end
