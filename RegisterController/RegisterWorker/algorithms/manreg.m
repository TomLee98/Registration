function status = manreg(movsrc_, movdst_, movtmpl_, regfrs_, regopt_)
%MANREG This function is two channels manual registration algorithm, which
%corresponds to the caller RegisterWorker
% Input:
%   - movsrc_: 1-by-1 regmov object, source registration movie
%   - movdst_: 1-by-1 regmov object, destination registration movie
%   - movtmpl_: m-by-n-by-p uint16 array, the registration template volume
%   - regfrs_: 1-by-n positive integer numeric array, indicating the frames
%              need to be aligned
%   - regopt_: 1-by-1 struct, with registration options
% Output:
%   - status: 1-by-1 double, 0 for normal exit, other code for bad exit
%
%   see also: regmov, regtmpl, fitgeotrans, imref2d, imwarp, cpselect

arguments
    movsrc_ (1,1)   regmov
    movdst_ (1,1)   regmov
    movtmpl_(:,:,:) uint16
    regfrs_ (1,1)   double {mustBePositive, mustBeInteger}
    regopt_ (1,1)   struct {mustBeRegistrationOption}
end

fmode = ["none", string(regfrs_)];

% extract variables because unix matlab bad memory management and memory leak
% new variables decrease total memory allocating times to one
dview = regopt_.DView;
proj = regopt_.Projection;
tf_type = regopt_.TformType;
interp_alg = regopt_.Interp;
tfenb = regopt_.TFEnable;
tfmat = regopt_.TFMatrix;
sc = regopt_.SC;
fc = regopt_.FC;
rs = regopt_.Resampling;
degree = regopt_.Degree;

% extract functional and structured channel data
% note that: avol with  dimension order X,Y,Z
avol_sc = grv(movsrc_, fmode, sc);
avol_fc = grv(movsrc_, fmode, fc);
pfunc = str2func(proj);

% extract filled value
fival_sc = mean(avol_sc(:,[1,end],:),"all");
fival_fc =  mean(avol_fc(:,[1,end],:),"all");

if tfenb == true
    % input matrix with higher priority
    rref = imref3d(size(avol_sc));
    if isMATLABReleaseOlderThan("R2022b")
        tf = affine3d(tfmat);
    else
        tf = affinetform3d(tfmat);
    end
    avol_sc = imwarp(avol_sc, tf, interp_alg, "OutputView",rref, ...
        "FillValues", fival_sc);
    avol_fc = imwarp(avol_fc, tf, interp_alg, "OutputView",rref, ...
        "FillValues", fival_fc);
else
    % get the projection on view direction
    switch dview
        case "XY"
            pjdim = 3;
        case "ZX"
            pjdim = 1;
        case "YZ"
            pjdim = 2;
        otherwise
    end

    switch lower(proj)
        case {'min', 'max'}
            pavol_sc = squeeze(pfunc(avol_sc, [], pjdim));
            ptmvol = squeeze(pfunc(movtmpl_, [], pjdim));
        case {'mean', 'median'}
            pavol_sc = squeeze(pfunc(avol_sc, pjdim));
            ptmvol = squeeze(pfunc(movtmpl_, pjdim));
        otherwise
    end

    % remove the possible over brighten pixels
    pamax = prctile(pavol_sc, 99.9,"all");
    pavol_sc(pavol_sc>pamax) = pamax;
    pamax = prctile(ptmvol, 99.9,"all");
    ptmvol(ptmvol>pamax) = pamax;

    % rescale the image to uint8
    pavol_sc = Remap(pavol_sc);
    ptmvol = Remap(ptmvol);

    % resampling the image for better alignment
    rs = rs(setdiff([2,1,3], pjdim, "stable"));
    pavol_sc_rz = imresize(pavol_sc, "Scale", rs, "Method","bilinear");
    ptmvol_rz = imresize(ptmvol, "Scale", rs, "Method","bilinear");

    % using cpselect to select the control points in GUI
    [mp, fp] = cpselect(pavol_sc_rz, ptmvol_rz, 'Wait',true);

    if isempty(mp) || isempty(fp)
        % manual exit
        status = 0;
        return;
    end

    % use cpcorr adjust the moving points
    mp = cpcorr(mp, fp, pavol_sc_rz, ptmvol_rz);

    if isMATLABReleaseOlderThan("R2022b")
        switch tf_type
            case "pwl"
                if size(fp, 1) < 4
                    tf = fewControlPoints(false);
                else
                    tf = fitgeotrans(mp, fp, "pwl");    % PiecewiseLinearTransformation2D
                end
            case "poly"
                if size(fp, 1) < nchoosek(degree+2, 2)
                    tf = fewControlPoints(false);
                else
                    tf = fitgeotrans(mp, fp, "polynomial", degree); % PolynomialTransformation2D
                end
            otherwise
                if tf_type == "affine"
                    if size(fp, 1) < 3
                        tf = fewControlPoints(false);
                    else
                        tf = fitgeotrans(mp, fp, "affine");    % affine2d
                    end
                else
                    if size(fp, 1) < 2
                        [mp, fp] = compensateTriangle(mp, fp, ptmvol_rz);
                    end

                    tf = fitgeotrans(mp, fp, "nonreflectivesimilarity");    % affine2d

                    if tf_type == "translation"
                        % omit scaling and rotation
                        tf.T(1:2,1:2) = eye(2);
                    elseif tf_type == "rigid"
                        % omit scaling
                        tf.T(1:2,1:2) = ...
                            tf.T(1:2,1:2) / sqrt(abs(det(tf.T(1:2,1:2))));
                    else
                        throw(MException("manreg:invalidGeometricTransformation", ...
                            "Invalid geometric transformation."));
                    end

                    % modify the translation because of resampling
                    tf.T(3,1:2) = tf.T(3,1:2)./fliplr(rs);
                end
        end
    else
        switch tf_type
            case "pwl"
                if size(fp, 1) < 4
                    tf = fewControlPoints(true);
                else
                    tf = fitgeotform2d(mp, fp, "pwl");      % PiecewiseLinearTransformation2D
                end
            case "poly"
                if size(fp, 1) < nchoosek(degree+2, 2)
                    tf = fewControlPoints(true);
                else
                    tf = fitgeotform2d(mp, fp, "polynomial", degree);   % 	PolynomialTransformation2D
                end
            otherwise
                if tf_type == "affine"
                    if size(fp, 1) < 3
                        tf = fewControlPoints(true);
                    else
                        tf = fitgeotform2d(mp, fp, "affine");    % affinetform2d
                    end
                    tf.A(1:2, 3) = tf.A(1:2, 3)./fliplr(rs)';
                else
                    if size(fp, 1) < 2
                        [mp, fp] = compensateTriangle(mp, fp, ptmvol_rz);
                    end

                    tf = fitgeotform2d(mp, fp, "similarity");    % simtform2d

                    if tf_type == "translation"
                        % omit scaling and rotation
                        tf.Scale = 1;
                        tf.RotationAngle = 0;
                    elseif tf_type == "rigid"
                        % omit scaling
                        tf.Scale = 1;
                    else
                        throw(MException("manreg:invalidGeometricTransformation", ...
                            "Invalid geometric transformation."));
                    end

                    % modify the translation because of resampling
                    tf.Translation = tf.Translation./fliplr(rs);
                end
        end
    end

    % create the output view
    rref = imref2d(size(ptmvol));

    % permute dimension for applying transformation autometically
    porder = [setdiff(1:3, pjdim, "stable"), pjdim];
    avol_sc = permute(avol_sc, porder);
    avol_fc = permute(avol_fc, porder);
    avol_sc = imwarp(avol_sc, tf, interp_alg, "OutputView",rref, ...
        "FillValues", fival_sc);
    avol_fc = imwarp(avol_fc, tf, interp_alg, "OutputView",rref, ...
        "FillValues", fival_fc);
    [~, porder] = ismember(1:3, porder);
    avol_sc = permute(avol_sc, porder);
    avol_fc = permute(avol_fc, porder);
end

% set the movdst
srv(movdst_, avol_sc, fmode, sc);
srv(movdst_, avol_fc, fmode, fc);
movdst_.Transformation{regfrs_, 3} = tf;

status = 0;
end

function [mp, fp] = compensateTriangle(mp, fp, img)
while true
    % add random points for constructing a triangle
    px = randi(size(img, 2), 3-size(fp,1), 1);
    py = randi(size(img, 1), 3-size(fp,1), 1);
    fpp = [fp; [px, py]];
    % avoid three points at a common line
    if abs(det(fpp(2:3,:)-fpp(1,:))) > eps
        mpp = [mp; mp+[px, py]-fp];
        fp = fpp;
        mp = mpp;
        break;
    end
end
end

function tf = fewControlPoints(new_flag)
if new_flag == false
    tf = affine2d();
else
    tf = affinetform2d();
end

warning("manreg:tooFewControlPointsPairs", ...
    "The number of control points does not meet the minimum requirements.");
end


function mustBeRegistrationOption(A)
if ~ismember("Mode", fieldnames(A))
    throw(MException("mustBeRegistrationOption:invalidOption", ...
        "Unsupported arguments input."));
end

switch A.Mode
    case "global"
        VALID_FIELD = ["Mode", "SubAlgorithm", "TformType", "Degree", "DView", "Interp", ...
            "Projection", "Resampling", "Isometric", "TFEnable", "TFMatrix", ...
            "SC", "FC", "Hardware"];
    otherwise
end

if ~all(ismember(fieldnames(A), VALID_FIELD))
    throw(MException("mustBeRegistrationOption:invalidOption", ...
        "Unsupported arguments input."));
end

end