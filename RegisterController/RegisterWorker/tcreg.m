function status = tcreg(movsrc_, movdst_, movtmpl_, regfrs_, regopt_)
%TCREG This function is two channels registration algorithm, which
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
%   see also: regmov, regtmpl, imregtform, imregcorr, imwarp, imregdemons, 
%   imhistmatchn, imref3d

arguments
    movsrc_ (1,1)   regmov
    movdst_ (1,1)   regmov
    movtmpl_(:,:,:) uint16
    regfrs_ (1,:)   double {mustBePositive, mustBeInteger}
    regopt_ (1,1)   struct {mustBeRegistrationOption}
end

% parallel checking
if isempty(gcp("nocreate"))
    throw(MException("tcreg:invalidParpoolStatus", ...
        "No workers are distributed."));
end

switch regopt_.Mode
    case "global"
        status = tcreg_global_cpu(movsrc_, movdst_, movtmpl_, regfrs_, regopt_);
    case "local"
        switch regopt_.Hardware
            case "cpu"
                status = tcreg_local_cpu(movsrc_, movdst_, movtmpl_, regfrs_, regopt_);
            case "cpu|gpu"
                status = tcreg_local_gpu(movsrc_, movdst_, movtmpl_, regfrs_, regopt_);
            otherwise
        end
    otherwise
        % just return
        status = 0;
end

end

function status = tcreg_global_cpu(movsrc, movdst, refvol, regfrs, regopt)
% This function use fast robust coarse registration pipeline
% for global motion estimation
% [Image Processing Pipeline]
%   (a)  downsampling
%   (b)  median filter (remove abnormal hardware noise) 
%        -> imhistmatchn (correct the global light change)
%   (c)  ->gaussian low pass filter(smooth for more robust)
%   (d)  imregmoment for initial estimation (>=2022b)
%   (d') imregcorr + z-direction optimization (<2022b)
%   (e)  upsampling
%   (f)  imregtform fine tuning

UPSAMPLE_FACTOR = 2;
TF = cell(numel(regfrs), 1);
fmode = ["none", string(regfrs).join(",")];
if regopt.DS == "auto"
    regds = inf;
else
    regds = 1/str2double(regopt.DS.extractBefore("X"));
end

% 1. extract the reference volume
[refvol_ds, ds_scale] = ReSample(refvol, regds);
[refvol_us, ~] = ReSample(refvol, UPSAMPLE_FACTOR);

% 2. extract functional and structured channel data
avol_sc = grv(movsrc, fmode, regopt.SC);
avol_fc = grv(movsrc, fmode, regopt.FC);

% 3. generate reference object and downsampling reference object
res_ds = [movsrc.MetaData.xRes/ds_scale, ...
          movsrc.MetaData.yRes/ds_scale, ...
          movsrc.MetaData.zRes/1];
res_us = [movsrc.MetaData.xRes/UPSAMPLE_FACTOR, ...
          movsrc.MetaData.yRes/UPSAMPLE_FACTOR, ...
          movsrc.MetaData.zRes/1];
rref = imref3d(size(refvol), movsrc.MetaData.xRes, movsrc.MetaData.yRes, movsrc.MetaData.zRes);
rref_ds  = imref3d(size(refvol_ds), res_ds(1), res_ds(2), res_ds(3));
rref_us = imref3d(size(refvol_us), res_us(1), res_us(2), res_us(3));

% extract local variables instead of struct spread
reg_modal = regopt.RegModal;
tf_type = regopt.TformType;
max_step = regopt.MaxStep;
min_step = regopt.MinStep;
iter_coeff = regopt.IterCoeff;
max_iter_n = regopt.MaxIterN;
vpl = regopt.VPL;
itpalg = regopt.Interp;

parfor m = 1:numel(regfrs)
    % downsampling  on selected volume
    avol_sc_m = avol_sc(:,:,:,m);
    avol_fc_m = avol_fc(:,:,:,m);
    avol_sc_m_ds = ReSample(avol_sc_m, ds_scale);
    avol_sc_m_ds = doPreProcessOn(avol_sc_m_ds, refvol_ds);
   
    % gaussian filter for robust motion correction
    avol_sc_m_ds = imgaussfilt3(avol_sc_m_ds, "FilterSize", [3,3,1], ...
        "Padding", "replicate");

    if isMATLABReleaseOlderThan("R2022b")
        % older than R2022b, call function estimate...
        % ptf: pre-transformation as affine3d object
        [ptf, ~] = imregopzr(avol_sc_m_ds, refvol_ds, res_ds);
        fival_sc = mean(avol_sc_m(:,[1,end],:),"all");
    else
        % call imregmoment for estimation
        % ptf: pre-transformation as affinetform3d object
        [ptf, ~] = imregmoment(avol_sc_m_ds, rref_ds, refvol_ds, rref_ds, ...
            "MedianThresholdBitmap",true);  % medianthreshold for more robust
        fival_sc = mean(avol_sc_m(:,[1,end],:),"all");
    end
    fival_fc =  mean(avol_fc_m(:,[1,end],:),"all");

    % 2 fold upsampling on selected volume
    avol_sc_m = doPreProcessOn(avol_sc_m, refvol);
    avol_sc_m_us = ReSample(avol_sc_m, 2);

    % set up the registration optimizer options
    [optimizer, metric] = imregconfig(reg_modal);
    switch reg_modal
        case "multimodal"
            optimizer.InitialRadius = max_step;
            optimizer.Epsilon = min_step;
            optimizer.GrowthFactor = iter_coeff;
        case "monomodal"
            optimizer.MaximumStepLength = max_step;
            optimizer.MinimumStepLength = min_step;
            optimizer.RelaxationFactor = iter_coeff;
        otherwise
    end
    optimizer.MaximumIterations = max_iter_n;

    % do fine registration: imregtform
    tf = imregtform(avol_sc_m_us, rref_us, refvol_us, rref_us, ...
        tf_type, optimizer, metric,...
        "PyramidLevels",vpl,...
        "InitialTransformation",ptf);

    % imwarp the transformation
    % replace avol_sc_m
    avol_sc_m = imwarp(avol_sc(:,:,:,m), rref, tf, ...
        itpalg,'OutputView',rref,'FillValues',fival_sc);
    avol_fc_m = imwarp(avol_fc(:,:,:,m), rref, tf, ...
        itpalg,'OutputView',rref,'FillValues',fival_fc);

    avol_sc(:,:,:,m) = avol_sc_m;
    avol_fc(:,:,:,m) = avol_fc_m;
    TF{m} = tf;
end

% set the movdst
srv(movdst, avol_sc, fmode, regopt.SC);
srv(movdst, avol_fc, fmode, regopt.FC);
movdst.Transformation(regfrs, 1) = TF;

status = 0;
end

function tcreg_local_cpu(movsrc, movdst, refvol, regfrs, regopt)

end

function tcreg_local_gpu(movsrc, movdst, refvol, regfrs, regopt)

end


% ====================== Utility Functions =======================
function vs = doPreProcessOn(vs_, refv_, mfsize_, hn_)
arguments
    vs_ 
    refv_   (:,:,:) double
    mfsize_ (1,3)   double {mustBeInteger, mustBePositive} = [3,3,3]
    hn_     (1,1)   double {mustBeMember(hn_, [256,512,1024,2048])} = 1024
end

vs_ = medfilt3(vs_, mfsize_, "replicate");
vs = imhistmatchn(vs_, refv_, hn_);
end

function [tf_est, movol_est] = imregopzr(movol_, refvol_, res_)
% This function use imregcorr and z optimization for transformation
% estimation on platform version < R2022b
arguments
    movol_  (:,:,:) double
    refvol_ (:,:,:) double
    res_    (1,3)   double      % [x,y,z] coordinate resolution, unit as um/pix
end

% do maximum z projection  for imregcorr
mov_img = max(movol_, [], 3);
ref_img = max(refvol_, [], 3);
rref = imref2d(size(ref_img), res_(1), res_(2));
rref3d = imref3d(size(refvol_), res_(1), res_(2), res_(3));
zlim = size(refvol_, 3)*[-1, 1];

% get imwarp filled value
fi_val = mean([mov_img(1,:)'; mov_img(end,:)'; mov_img(:,1); mov_img(:,end)]);

% 2-D rigid transformation estimation
tf2d_ = imregcorr(mov_img, rref, ref_img, rref, "rigid");  % rigid2d object

% optimize the z shift by immse as loss function
optf = @(x)opfun(x, movol_, refvol_, rref3d, tf2d_, fi_val);
fminbnd_opts = optimset('MaxFunEvals',500, 'MaxIter',100, 'TolX', 1e-3);
[z_, fval] = fminbnd(optf, zlim(1), zlim(2), fminbnd_opts);

% transform rigid2d object to affine3d object as imregtform initialized
% transformation estimation
TMat = [[[tf2d_.T(1:2,1:2),[0;0]];[0,0,1]];[tf2d_.T(3,1:2),z_]];
TMat = [TMat,[0;0;0;1]];
tf_est = affine3d(TMat);

if nargout == 3
    % 1 memory copy from <imwarp>
    movol_est = imwarp(movol_, rref3d, tf_est, "linear",...
        "OutputView",rref3d, 'FillValues',fi_val);
end

    function f = opfun(z_, mov_, ref_, ra_, tf_, fival_)
        rot = [[tf_.Rotation;[0,0]], [0,0,1]];  % no z correlated rotation
        trans = [tf_.Translation, z_];          % add shifts on z dimension 
        T = rigid3d(rot, trans);
        % imrotate3 and imtranslate?
        mov_ = imwarp(mov_, ra_, T, "linear", ...
            "FillValues",fival_, "OutputView",ra_);
        f = immse(mov_, ref_);
    end

end

function mustBeRegistrationOption(A)
VALID_FIELD = ["Mode", "RegModal", "TformType", "MaxStep", "MinStep", "MaxIterN", ...
    "IterCoeff", "VPL", "Interp", "CoRegC", "DS", "SC", "FC", "Hardware"];

if ~all(ismember(fieldnames(A), VALID_FIELD))
    throw(MException("mustBeRegistrationOption:invalidOption", ...
        "Unsupported arguments input."));
end

% Escape value checking for tcreg calling faster

end