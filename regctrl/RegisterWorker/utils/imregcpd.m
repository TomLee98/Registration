function [dfield, regVol] = imregcpd(moving, fixed, options)
%IMREGCPD This function use point cloud algorithm for fast image non-rigid
%registration

% Copyright 2025 Weihan, Li.

% Validate the input images.

arguments
    moving {mustBeNumeric, mustBeNonNan, mustBeFinite,...
        mustBeNonsparse, mustBeNonempty, mustBeReal}
    fixed (:,:,:) {mustBeNumeric, mustBeNonNan, mustBeFinite,...
        mustBeNonsparse, mustBeNonempty, mustBeReal, mustHaveSameSize(moving, fixed)}

    options.DenoiseThreshold (1,1) {mustBeNumeric, mustBeNonsparse, ...
        mustBeNonempty, mustBeReal, mustBePositive} = 1;

    options.GridSpacing (1,1) {mustBeNumeric, mustBeNonsparse, mustBeFinite, ...
        mustBeNonempty, mustBeReal, mustBePositive} = 3;

    options.GridRegularization (1,1) {mustBeNumeric, mustBeNonsparse, mustBeFinite, ...
        mustBeNonempty, mustBeReal, mustBePositive, mustBeInRange(options.GridRegularization, 1.5, 3)} = 2;
    options.OutlierRatio (1,1) {mustBeNumeric, mustBeNonsparse, mustBeFinite, ...
        mustBeNonempty, mustBeReal, mustBeInRange(options.OutlierRatio, 0, 1)} = 0.1;
    options.MaxIterations (1,1) {mustBeNumeric, mustBeNonsparse, mustBeFinite, ...
        mustBeNonempty, mustBeReal, mustBePositive, mustBeInteger} = 20;
    options.Verbose (1,1) {mustBeNumericOrLogical, mustBeNonsparse, mustBeFinite, ...
        mustBeNonempty, mustBeReal} = false;

    options.InterpMethod (1,1) {mustBeNonzeroLengthText, ...
        mustBeMember(options.InterpMethod, ["linear","nearest","natural"])} = "linear";
end

grid_step = options.GridSpacing;
th = options.DenoiseThreshold;
outlier_ratio = options.OutlierRatio;
max_iter_n = options.MaxIterations;
grid_regular = options.GridRegularization;
verbose = options.Verbose;
itp_method = options.InterpMethod;

% split as foreground and background
vol_bw = imsegkmeans3(moving, 2); vol_bw = logical(vol_bw - 1);
ref_bw = imsegkmeans3(fixed, 2);  ref_bw = logical(ref_bw - 1);

% translate to point cloud and apply cpd algorithm
[py, px, pz] = ind2sub(size(ref_bw), find(ref_bw));
ref_pt = pointCloud([px, py, pz]);
[py, px, pz] = ind2sub(size(vol_bw), find(vol_bw));
vol_pt = pointCloud([px, py, pz]);

% remove outliers (may be slow, optional)
if ~isinf(th)
    ref_pt = pcdenoise(ref_pt, "Threshold", th);
    vol_pt = pcdenoise(vol_pt, "Threshold", th);
end

% downsampling for fast and robust registration
ref_pt = pcdownsample(ref_pt, "gridAverage", grid_step);    % for pixel binning size
vol_pt = pcdownsample(vol_pt, "gridAverage", grid_step);

% estimate the non-rigid transformation
tf = pcregistercpd(vol_pt, ref_pt, "Transform","Nonrigid", ...
    "OutlierRatio",outlier_ratio, "MaxIterations", max_iter_n, ...
    "InteractionSigma",grid_regular, "Verbose", verbose);

% griddata and apply to points convex hull
[X, Y, Z] = meshgrid(0.5:size(ref_bw,2)-0.5, 0.5:size(ref_bw,1)-0.5, 0.5:size(ref_bw,3)-0.5);
dFx = griddata(vol_pt.Location(:,1), vol_pt.Location(:,2), vol_pt.Location(:,3), ...
    tf(:,1), X, Y, Z, itp_method);
dFx(isnan(dFx)) = 0;
dFy = griddata(vol_pt.Location(:,1), vol_pt.Location(:,2), vol_pt.Location(:,3), ...
    tf(:,2), X, Y, Z, itp_method);
dFy(isnan(dFy)) = 0;
dFz = griddata(vol_pt.Location(:,1), vol_pt.Location(:,2), vol_pt.Location(:,3), ...
    tf(:,3), X, Y, Z, itp_method);
dFz(isnan(dFz)) = 0;

% concatnate for displacement field
dfield = cat(4, dFx, dFy, dFz);

if nargout == 2
    % apply imwarp
    regVol = imwarp(moving, dfield);
end

end

function mustHaveSameSize(moving, fixed)

[rf, cf, chf] = size(fixed);
validateattributes(moving, {'numeric'}, {'size', [rf cf chf]});
end