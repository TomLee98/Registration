function status = tcreg(movsrc_, movdst_, movtmpl_, regfrs_, regopt_)
%TCREG This function is two channels registration algorithm, which
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
%   see also: regmov, regtmpl, regopt, imregcoarse, imregtform, imwarp, 
%             imregdemons, imhistmatchn, imref3d

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

% ======================== GLOBAL ALGORITHM ===========================
% This function use fast robust coarse registration pipeline
% for global motion estimation
function status = tcreg_global_cpu(movsrc, movdst, refvol, regfrs, regopt)

TF = cell(numel(regfrs), 1);
fmode = ["none", string(regfrs).join(",")];
if regopt.DS == "auto"
    regds = inf;
else
    regds = 1/str2double(regopt.DS.extractBefore("X"));
end

mf = regopt.MedianFilter;
df = regopt.DilateFilter;
dfh = regopt.DilateFilterEnh;
gf = regopt.GaussianFilter;
ga = regopt.Gamma;

% 1. extract the reference volume
[ds_scale, refvol_ds] = Resample(refvol, regds);
if dfh == true
    [refvol_ds_pp, roiref] = preproc_tc(refvol_ds, df, mf, gf, ga);
    if ~isempty(roiref)
        roiref(:, 1:3) = [];    % keep w,h,d
    else
        if "mmt"==regopt.CoarseAlg
            throw(MException("tcreg:invalidDilateFilterArgs", ...
                "No valid object in volume under the given arguments."));
        end
    end
else
    roiref = double.empty(0, 3);
    refvol_ds_pp = preproc_tc(refvol_ds, df, mf, gf, ga);
end

% 2. extract functional and structured channel data
avol_sc = grv(movsrc, fmode, regopt.SC);
avol_fc = grv(movsrc, fmode, regopt.FC);

% 3. generate reference object and downsampling reference object
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
max_shift_z = regopt.MaxZOptShift;
zopt_tol = regopt.TolZOpt;
vpl = regopt.glVPL;
itpalg = regopt.Interp;

parfor m = 1:numel(regfrs)
    % down sampling on selected volume
    avol_sc_m = avol_sc(:,:,:,m);
    avol_fc_m = avol_fc(:,:,:,m);
    [~, avol_sc_m_ds] = Resample(avol_sc_m, ds_scale);

    if dfh == true
        avol_sc_m_ds_pp = preproc_tc(avol_sc_m_ds, df, mf, gf, ga, roiref);
    else
        avol_sc_m_ds_pp = preproc_tc(avol_sc_m_ds, df, mf, gf, ga);
    end

    % use imregcoarse for better initialized transformation
    % where the preprocess volumes are needed
    ptf = imregcoarse(avol_sc_m_ds_pp, refvol_ds_pp, res_ds, true, ...
            max_shift_z, zopt_tol, coarse_alg, coarse_args);

    fival_sc = mean(avol_sc_m(:,[1,end],:),"all");
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
    tf = imregtform(avol_sc_m_ds, rref_ds, refvol_ds, rref_ds, tf_type, ...
        optimizer, metric, "PyramidLevels",vpl, "InitialTransformation",ptf);

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
%======================================================================

%======================= LOCAL ALGORITHM (CPU)=========================
% This function use fast robust coarse registration pipeline
% for local motion estimation by using cpu
function status = tcreg_local_cpu(movsrc, movdst, refvol, regfrs, regopt)

DF = cell(numel(regfrs), 1);
fmode = ["none", string(regfrs).join(",")];

% 1. extract the brightness reference volume
[~, refvol_ds] = Resample(refvol, 1/4);

% 2. extract functional and structured channel data
avol_sc = grv(movsrc, fmode, regopt.SC);
avol_fc = grv(movsrc, fmode, regopt.FC);

% extract local variables instead of struct spread
pr = [movsrc.MetaData.xRes, movsrc.MetaData.yRes, movsrc.MetaData.zRes];
subalg = regopt.SubAlgorithm;
max_itern = regopt.MaxIterN;
afs = regopt.AFS;
gr = regopt.GR;
gs = regopt.GS;
if subalg == "usual"
    vpl = regopt.dfVPL;
elseif subalg == "advanced"
    vpl = regopt.dmVPL;
end
itpalg = regopt.Interp;
img_rehist = regopt.ImageRehist;
repacc = regopt.RepAcc;

parfor m = 1:numel(regfrs)
    % down sampling  on selected volume
    avol_sc_m = avol_sc(:,:,:,m);
    avol_fc_m = avol_fc(:,:,:,m);

    if img_rehist == true
        avol_sc_m = imhistmatchn(avol_sc_m, refvol_ds, repacc);
    end

    if subalg == "usual"
        % parameters: GR, GL, VPL are useful
        [df, ~] = imregdeform(avol_sc_m, refvol, "GridSpacing",gs, ...
            "GridRegularization",gr,"NumPyramidLevels",vpl, ...
            "PixelResolution",pr,"DisplayProgress",false);
    elseif subalg == "advanced"
        [df, ~] = imregdemons(avol_fc_m, refvol,...
        max_itern, "AccumulatedFieldSmoothing", afs,...
        "PyramidLevels", vpl, "DisplayWaitbar",false);
    else
        df = []; %#ok<NASGU>
        throw(MException("tcreg:local:invalidSubAlgorithm", ...
            "Unregistered sub algorithm."));
    end
    fival_sc = mean(avol_sc_m(:,[1,end],:),"all");
    fival_fc =  mean(avol_fc_m(:,[1,end],:),"all");

    % imwarp the displacement field
    % replace avol_sc_m
    avol_sc_m = imwarp(avol_sc(:,:,:,m), df, itpalg, "FillValues",fival_sc);
    avol_fc_m = imwarp(avol_fc(:,:,:,m), df, itpalg, "FillValues",fival_fc);

    avol_sc(:,:,:,m) = avol_sc_m;
    avol_fc(:,:,:,m) = avol_fc_m;
    DF{m} = df;
end

% set the movdst
srv(movdst, avol_sc, fmode, regopt.SC);
srv(movdst, avol_fc, fmode, regopt.FC);
movdst.Transformation(regfrs, 2) = DF;

status = 0;
end
%=======================================================================

%======================= LOCAL ALGORITHM (GPU)=========================
% This function use fast robust coarse registration pipeline
% for global motion estimation by using gpu
function status = tcreg_local_gpu(movsrc, movdst, refvol, regfrs, regopt)
DF = cell(numel(regfrs), 1);
fmode = ["none", string(regfrs).join(",")];

% 1. extract the brightness reference volume
[~, refvol_ds] = Resample(refvol, 1/4);

% 2. extract functional and structured channel data
avol_sc = grv(movsrc, fmode, regopt.SC);
avol_fc = grv(movsrc, fmode, regopt.FC);

% extract local variables instead of struct spread
pr = [movsrc.MetaData.xRes, movsrc.MetaData.yRes, movsrc.MetaData.zRes];
subalg = regopt.SubAlgorithm;
max_itern = regopt.MaxIterN;
afs = regopt.AFS;
gr = regopt.GR;
gs = regopt.GS;
vpl = regopt.VPL;
itpalg = regopt.Interp;
img_rehist = regopt.ImageRehist;
repacc = regopt.RepAcc;

parfor m = 1:numel(regfrs)
    % down sampling  on selected volume
    avol_sc_m = avol_sc(:,:,:,m);
    avol_fc_m = avol_fc(:,:,:,m);

    if img_rehist == true
        avol_sc_m = imhistmatchn(avol_sc_m, refvol_ds, repacc);
    end

    if subalg == "usual"
        % imregdeform only supports cpu calculation
        % parameters: GR, GL, VPL are useful
        [df, ~] = imregdeform(avol_sc_m, refvol, "GridSpacing",gs, ...
            "GridRegularization",gr,"NumPyramidLevels",vpl, ...
            "PixelResolution",pr,"DisplayProgress",false);
    elseif subalg == "advanced"
        % spread memory variable to GPU memory
        refvol_ga = gpuArray(refvol);
        avol_sc_m_ga = gpuArray(avol_sc_m);

        [df_ga, ~] = imregdemons(avol_sc_m_ga, refvol_ga,...
        max_itern, "AccumulatedFieldSmoothing", afs,...
        "PyramidLevels", vpl, "DisplayWaitbar",false);
        df = gather(df_ga);
    else
        df = []; %#ok<NASGU>
        throw(MException("tcreg:localgpu:invalidSubAlgorithm", ...
            "Unregistered sub algorithm."));
    end
    fival_sc = mean(avol_sc_m(:,[1,end],:),"all");
    fival_fc =  mean(avol_fc_m(:,[1,end],:),"all");

    % imwarp the displacement field
    % replace avol_sc_m, cpu running
    avol_sc_m = imwarp(avol_sc(:,:,:,m), df, itpalg, "FillValues",fival_sc);
    avol_fc_m = imwarp(avol_fc(:,:,:,m), df, itpalg, "FillValues",fival_fc);

    avol_sc(:,:,:,m) = avol_sc_m;
    avol_fc(:,:,:,m) = avol_fc_m;
    DF{m} = df;
end

% set the movdst
srv(movdst, avol_sc, fmode, regopt.SC);
srv(movdst, avol_fc, fmode, regopt.FC);
movdst.Transformation(regfrs, 2) = DF;

status = 0;
end


% ====================== Utility Functions =======================

function mustBeRegistrationOption(A)
if ~ismember("Mode", fieldnames(A))
    throw(MException("mustBeRegistrationOption:invalidOption", ...
        "Unsupported arguments input."));
end

VALID_FIELD_PUBLIC = ["Mode", "SubAlgorithm", "SC", "FC", "Hardware", ...
    "MaxIterN", "Interp"];

switch A.Mode
    case "global"
        VALID_FIELD_PRIVATE = ["RegModal", "AreaMask", "MedianFilter", "GaussianFilter", ...
            "DilateFilter", "MaxZOptShift", "TolZOpt", "Gamma", "TformType", ...
            "MaxStep", "MinStep",  "IterCoeff", "DS", "CoarseAlg", "CoarseArgs",...
            "DilateFilterEnh", "glVPL"];
    case "local"
        VALID_FIELD_PRIVATE = ["dmVPL", "dfVPL", "AFS", "GR", "GS", "ImageRehist", "RepAcc"];
    otherwise
end

if ~all(ismember(fieldnames(A), [VALID_FIELD_PUBLIC, VALID_FIELD_PRIVATE]))
    throw(MException("mustBeRegistrationOption:invalidOption", ...
        "Unsupported arguments input."));
end

% Escape value checking for tcreg running faster

end