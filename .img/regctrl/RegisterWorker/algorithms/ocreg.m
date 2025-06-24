function status = ocreg(movsrc_, movdst_, movtmpl_, regfrs_, regopt_)
%OCREG This function is one channels registration algorithm, which
%corresponds to the caller RegisterWorker
% Input:
%   - movsrc_: 1-by-1 regmov object, source registration movie
%   - movdst_: 1-by-1 regmov object, destination registration movie
%   - movtmpl_:1-by-1 uint16 array, the registration template volume
%   - regfrs_: 1-by-n positive integer numeric array, indicating the frames
%              need to be aligned
%   - regopt_: 1-by-1 struct, with registration options
% Output:
%   - status: 1-by-1 double, 0 for normal exit, other code for bad exit
%
%   see also: regmov, regtmpl, imregtform, imregcorr, imwarp, imref3d

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
        status = ocreg_global_cpu(movsrc_, movdst_, movtmpl_, regfrs_, regopt_);
    otherwise
        % just return
        status = 0;
end

end

% ======================== GLOBAL ALGORITHM ===========================
% This function use fast robust coarse registration pipeline
% for global motion estimation
function status = ocreg_global_cpu(movsrc, movdst, refvol, regfrs, regopt)
TF = cell(numel(regfrs), 1);
fmode = ["none", string(regfrs).join(",")];
if regopt.DS == "auto"
    regds = inf;
else
    regds = 1/str2double(regopt.DS.extractBefore("X"));
end

mf = regopt.MedianFilter;
gf = regopt.GaussianFilter;

% 1. extract the reference volume
[ds_scale, refvol_ds] = Resample(refvol, regds);
refvol_ds_pp = preproc_oc(refvol_ds, mf, gf);

% 2. extract functional and structured channel data
avol_fc = grv(movsrc, fmode, regopt.FC);

% 3. generate reference object and down sampling reference object
res_ds = [movsrc.MetaData.xRes/ds_scale, ...
          movsrc.MetaData.yRes/ds_scale, ...
          movsrc.MetaData.zRes/1];
rref = imref3d(size(refvol), movsrc.MetaData.xRes, movsrc.MetaData.yRes, movsrc.MetaData.zRes);
rref_ds  = imref3d(size(refvol_ds), res_ds(1), res_ds(2), res_ds(3));

% extract local variables instead of struct spread
reg_modal = regopt.RegModal;
tf_type = regopt.TformType;
coarse_alg = regopt.CoarseAlg;
coarse_args = regopt.CoarseArgs;
max_step = regopt.MaxStep;
min_step = regopt.MinStep;
iter_coeff = regopt.IterCoeff;
max_itern = regopt.MaxIterN;
zopt_tol = regopt.TolZOpt;
vpl = regopt.VPL;
itpalg = regopt.Interp;

parfor m = 1:numel(regfrs)
    % down sampling on selected volume
    avol_fc_m = avol_fc(:,:,:,m);
    [~, avol_fc_m_ds] = Resample(avol_fc_m, ds_scale);

    avol_fc_m_ds_pp = preproc_oc(avol_fc_m_ds, mf, gf);

    % use imregcoarse for better initialized transformation
    % where the preprocess volumes are needed
    ptf = imregcoarse(avol_fc_m_ds_pp, refvol_ds_pp, res_ds, false, ...
            zopt_tol, coarse_alg, coarse_args);

    fival_fc =  mean(avol_fc_m(:,[1,end],:),"all");

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
    optimizer.MaximumIterations = max_itern;

    % do fine registration: imregtform base on intensity
    tf = imregtform(avol_fc_m_ds, rref_ds, refvol_ds, rref_ds, tf_type, ...
        optimizer, metric, "PyramidLevels",vpl, "InitialTransformation",ptf);

    % imwarp the transformation
    % replace avol_fc_m
    avol_fc_m = imwarp(avol_fc(:,:,:,m), rref, tf, ...
        itpalg,'OutputView',rref,'FillValues',fival_fc);

    avol_fc(:,:,:,m) = avol_fc_m;
    TF{m} = tf;
end

% set the movdst
srv(movdst, avol_fc, fmode, regopt.FC);
movdst.Transformation(regfrs, 1) = TF;

status = 0;
end

function mustBeRegistrationOption(A)
if ~ismember("Mode", fieldnames(A))
    throw(MException("mustBeRegistrationOption:invalidOption", ...
        "Unsupported arguments input."));
end

VALID_FIELD = ["Mode", "SubAlgorithm", "RegModal", "TformType", "GaussianFilter", ...
    "MedianFilter", "MaxStep", "MinStep", "MaxIterN", "CoarseAlg", "CoarseArgs", ...
    "TolZOpt", "IterCoeff", "glVPL", "Interp", "DS", "SC", "FC", "Hardware"];

if ~all(ismember(fieldnames(A), VALID_FIELD))
    throw(MException("mustBeRegistrationOption:invalidOption", ...
        "Unsupported arguments input."));
end

% Escape value checking for ocreg calling faster

end