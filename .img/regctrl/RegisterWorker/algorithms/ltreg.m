function status = ltreg(movsrc_, movdst_, movtmpl_, regfrs_, regopt_)
%LTREG This function is long term registration algorithm, which
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

TF = cell(numel(regfrs_), 1);
fmode = ["none", string(regfrs_).join(",")];
if regopt_.DS == "auto"
    regds = inf;
else
    regds = 1/str2double(regopt_.DS.extractBefore("X"));
end

% 1. extract the reference volume
[ds_scale, refvol_ds] = Resample(movtmpl_, regds);

% 2. extract functional and structured channel data
avol_sc = grv(movsrc_, fmode, regopt_.SC);
avol_fc = grv(movsrc_, fmode, regopt_.FC);

% 3. generate reference object and down sampling reference object
rs = [movsrc_.MetaData.xRes, movsrc_.MetaData.yRes, movsrc_.MetaData.zRes];
rs_ds = [movsrc_.MetaData.xRes/ds_scale, ...
          movsrc_.MetaData.yRes/ds_scale, ...
          movsrc_.MetaData.zRes/1];
rref = imref3d(size(movtmpl_), rs(1), rs(2), rs(3));
rref_ds  = imref3d(size(refvol_ds), rs_ds(1), rs_ds(2), rs_ds(3));

% extract local variables instead of struct spread
max_step = regopt_.MaxStep;
min_step = regopt_.MinStep;
iter_coeff = regopt_.IterCoeff;
max_itern = regopt_.MaxIterN;
dfsize = regopt_.DilateFilter;
itpalg = regopt_.Interp;
vpl = fix(log(movsrc_.MetaData.slices)/log(4)) + 1;
rc = regopt_.RegChain;  % registration chain object

[tfs, kvs, ids] = rc.vacquire(regfrs_);

parfor m = 1:numel(regfrs_)
    % set up the registration optimizer options
    [optimizer, metric] = imregconfig("monomodal");
    optimizer.MaximumStepLength = max_step;
    optimizer.MinimumStepLength = min_step;
    optimizer.RelaxationFactor = iter_coeff;
    optimizer.MaximumIterations = max_itern;

    % loading to parallel worker
    kvs_par = kvs{ids(m)}; %#ok<PFBNS>
    tfs_par = tfs{ids(m)}; %#ok<PFBNS>

    % STEP 1: down sampling and run localized registration to key frame
    % selected volume
    avol_sc_m = avol_sc(:,:,:,m);
    avol_fc_m = avol_fc(:,:,:,m);

    fival_sc = mean(avol_sc_m(:,[1,end],:),"all");
    fival_fc =  mean(avol_fc_m(:,[1,end],:),"all");

    % preprocess for robust registration
    avol_sc_proc_m = preproc_tc(avol_sc_m, dfsize);   % dilate filter for robust estimation
    kvs_par = preproc_tc(kvs_par, dfsize);

    [~, avol_sc_m_ds] = Resample(avol_sc_proc_m, ds_scale);
    [~, kvol_m_ds] = Resample(kvs_par, ds_scale);

    % use imregcoarse for better initialized transformation
    % where the preprocess volumes are needed
    ptf = imregcoarse(avol_sc_m_ds, kvol_m_ds, rs_ds, true);

    % do fine registration, align to key frame
    ptf = imregtform(avol_sc_m_ds, rref_ds, kvol_m_ds, rref_ds, "affine", ...
        optimizer, metric, "PyramidLevels",vpl, "InitialTransformation",ptf);

    % STEP 2: use pre-aligned chain information as initial transformation
    ptf.A = tfs_par*ptf.A;

    % refine the registration, align to template volume
    % use ptf.A as estimation because of world coordinate system enabled
    tf = imregtform(avol_sc_m, rref, movtmpl_, rref, "affine", ...
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

% set the movdst_
srv(movdst_, avol_sc, fmode, regopt_.SC);
srv(movdst_, avol_fc, fmode, regopt_.FC);
movdst_.Transformation(regfrs_, 1) = TF;

status = 0;
end


function mustBeRegistrationOption(A)
VALID_FIELD = ["Mode", "SubAlgorithm", "Keyframes", "AutoKeyframe", "AutoTemplate", ...
    "TGridMinMax", "RegChain", "MaxZOptShift", "DilateFilter", "DS", "Interp", ...
    "MaxStep", "MinStep", "MaxIterN", "IterCoeff", "SC", "FC", "Hardware"];

if ~all(ismember(fieldnames(A), VALID_FIELD))
    throw(MException("mustBeRegistrationOption:invalidOption", ...
        "Unsupported arguments input."));
end

if isempty(A.RegChain) || ~isvalid(A.RegChain)
    throw(MException("mustBeRegistrationOption:invalidOption", ...
        "Registration chain is empty or invalid."));
end

% Escape some value checking for tcreg calling faster

end