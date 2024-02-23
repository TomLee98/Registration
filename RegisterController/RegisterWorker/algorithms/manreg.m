function status = manreg(movsrc_, movdst_, movtmpl_, regfrs_, regopt_)
%MANREG This function is two channels manual registration algorithm, which
%corresponds to the caller RegisterWorker
% Input:
%   - movsrc_: 1-by-1 regmov object, source registration movie
%   - movdst_: 1-by-1 regmov object, destination registration movie
%   - movtmpl_: m-by-n-by-p uint16 array, the registration template volume
%   - regfrs_: 1-by-n positive integer numeric array, indicating the frames
%              need to be aligned
%   - regopt_: 1-by-1 regopt object, with registration options
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

% extract functional and structured channel data
% note that: avol with  dimension order X,Y,Z
avol_sc = grv(movsrc_, fmode, regopt_.SC);
avol_fc = grv(movsrc_, fmode, regopt_.FC);
pfunc = str2func(regopt_.Projection);

% get the projection on view direction
switch regopt_.DView
    case "XY"
        pjdim = 3;
    case "ZX"
        pjdim = 1;
    case "YZ"
        pjdim = 2;
    otherwise
end

switch lower(regopt_.Projection)
    case {'min', 'max'}
        pavol_sc = squeeze(pfunc(avol_sc, [], pjdim));
        ptmvol = squeeze(pfunc(movtmpl_, [], pjdim));
    case {'mean', 'median'}
        pavol_sc = squeeze(pfunc(avol_sc, pjdim));
        ptmvol = squeeze(pfunc(movtmpl_, pjdim));
    otherwise
end

% extract filled value
fival_sc = mean(avol_sc(:,[1,end],:),"all");
fival_fc =  mean(avol_fc(:,[1,end],:),"all");

% remove the possible over brighten pixels
pamax = prctile(pavol_sc, 99.9,"all");
pavol_sc(pavol_sc>pamax) = pamax;
pamax = prctile(ptmvol, 99.9,"all");
ptmvol(ptmvol>pamax) = pamax;

% rescale the image to uint8
pavol_sc = Remap(pavol_sc);
ptmvol = Remap(ptmvol);

% resampling the image for better alignment
rs = regopt_.Resampling(setdiff([2,1,3], pjdim, "stable"));
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
    switch regopt_.TformType
        case "pwl"
            tf = fitgeotrans(mp, fp, "pwl");
        case "poly"
            tf = fitgeotrans(mp, fp, "polynomial", regopt_.Degree);
        otherwise
            if size(fp,1) < 3
                while true
                    % add random points for constructing a triangle
                    px = randi(size(ptmvol_rz, 2), 3-size(fp,1), 1);
                    py = randi(size(ptmvol_rz, 1), 3-size(fp,1), 1);
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
            tf = fitgeotrans(mp, fp, "nonreflectivesimilarity");    % affine2d
            if regopt_.TformType == "translation"
                % omit scaling and rotation
                tf.T(1:2,1:2) = eye(2);     
            elseif regopt_.TformType == "rigid"
                % omit scaling
                tf.T(1:2,1:2) = ...
                    tf.T(1:2,1:2) / sqrt(abs(det(tf.T(1:2,1:2))));
            else
            end
            % modify the translation because of resampling
            tf.T(3,1:2) = tf.T(3,1:2)./fliplr(rs);
    end
else
    switch regopt_.TformType
        case "pwl"
            tf = fitgeotrans(mp, fp, "pwl");
        case "poly"
            tf = fitgeotrans(mp, fp, "polynomial", regopt_.Degree);
        otherwise
            if size(fp,1) < 3
                while true
                    % add random points for constructing a triangle
                    px = randi(size(ptmvol_rz, 2), 3-size(fp,1), 1);
                    py = randi(size(ptmvol_rz, 1), 3-size(fp,1), 1);
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
            tf = fitgeotform2d(mp, fp, "similarity");    % simtform2d
            if regopt_.TformType == "translation"
                % omit scaling and rotation
                tf.Scale = 1;
                tf.RotationAngle = 0;
            elseif regopt_.TformType == "rigid"
                % omit scaling
                tf.Scale = 1;
            else
            end

            % modify the translation because of resampling
            tf.Translation = tf.Translation./fliplr(rs);
    end
end

% create the output view
rref = imref2d(size(ptmvol));

% permute dimension for applying transformation autometically
porder = [setdiff(1:3, pjdim, "stable"), pjdim];
avol_sc = permute(avol_sc, porder);
avol_fc = permute(avol_fc, porder);
avol_sc = imwarp(avol_sc, tf, regopt_.Interp, "OutputView",rref, ...
    'FillValues', fival_sc);
avol_fc = imwarp(avol_fc, tf, regopt_.Interp, "OutputView",rref, ...
    'FillValues', fival_fc);
[~, porder] = ismember(1:3, porder);
avol_sc = permute(avol_sc, porder);
avol_fc = permute(avol_fc, porder);

% set the movdst
srv(movdst_, avol_sc, fmode, regopt_.SC);
srv(movdst_, avol_fc, fmode, regopt_.FC);
movdst_.Transformation{regfrs_, 3} = tf;

status = 0;
end


function mustBeRegistrationOption(A)
if ~ismember("Mode", fieldnames(A))
    throw(MException("mustBeRegistrationOption:invalidOption", ...
        "Unsupported arguments input."));
end

switch A.Mode
    case "global"
        VALID_FIELD = ["Mode", "SubAlgorithm", "TformType", "Degree", "DView", "Interp", ...
            "Projection", "Resampling", "Isometric", "SC", "FC", "Hardware"];
    otherwise
end

if ~all(ismember(fieldnames(A), VALID_FIELD))
    throw(MException("mustBeRegistrationOption:invalidOption", ...
        "Unsupported arguments input."));
end

% Escape value checking for tcreg calling faster

end