function [signal,tform,marker] = register3D_tc_auto(varargin)
%REGISTER - This function help to align 3D volume series
%   This function can automatically align two channels volumetric images
%
%   signal = register3D()
%   [signal,tform,marker] = register3D()
%   [signal,tform,marker] = register3D(filename)
%   [signal,tform,marker] = register3D(movinfo,filename)
%   [signal,tform,marker] = register3D(___,Name,Value)
%
%   input:
%   - (optional) movinfo:  the structure of movie data(4/5-D) and
%                          information structure with field {'mov','opts'}
%   - (optional) filename: the name of file need to be aligned, only *.ims,
%                          *.nd2 and *.tiff is valid, image dimension: x*y*(c)*z*t
%   - varargin:
%       - RefVol: the fixed volume indicator, struct with {'G','L'}, 'G' for
%                 'global' and 'L' for 'local', each element is 1*2 string,
%                 for [method, range]. method indicates Ref calculate method,
%                 range indicates the sampling frames.
%                 Example: RefVol = struct("G",["mean","1,2,3"], "L",["mean","1,2,3"]);
%       - RegFrames: the frames series need to be registered
%       - RegMode: the register mode, ["global-only"(default),"local-only","mix"]
%       - Modal: the modal of aligned pair, ["multimodal","monomodal"(default)]
%       - Tform: the transform method when aligned, ["translation"(default),"rigid",
%                   "similarity","affine"]
%       - AutoContrast: the flag for auto contrast on marker channel, for
%                       better local registration, true(default)/false
%       - CompAcc: the compensate accuracy, the imhistmatchn binning counts
%       - InitStep:  the initial step when optimizal two volumes loss, it is
%                   initial radius(default: 6.25e-3) when using "multimodal", or
%                   maximum step length(default: 6.25e-2) when using "monomodal"
%       - MinStep: the minimum iteration step when optimizal two volumes
%                  loss, which meature the convergance
%       - MaxIterN_Rigid: the maximum iteration number when alignment, "auto" (default)
%       - MaxIterN_NonRigid: the maximum iteration number when alignment, 100 (default)
%       - IterCoeff: the iteration coefficient for loop, if Modal is monomodal,
%                    then the range of IterCoeff is (0,1), 0.5 as default;
%                    and if Modal is multimodal, then the range of
%                    IterCoeff is (0,+inf), 1.05 as default
%       - AFS: the accumulative field smoothing, the AFS bigger and the
%              displacement is more smooth, usual 0.5-3, 1 as default
%       - CoRegC: the clearity of objects when running coarse registration,
%                 small number indicate robust but information loss, 0-100, 50(default)
%       - DS: the downsampling options for speed up and more robust, could
%         be "auto", "2X2", "3X3", "4X4", "auto" as default
%       - Chla: the channel you want to align, this section is enabled
%                only if multi-channel ["r"(default),"g","b"]
%       - Chls: the channel you want to extract signal ["r","g"(default),"b"]
%       - VPL_Rigid: the valid pyramid levels for imregtform, 3 (default)
%       - VPL_NonRigid: the valid pyramid levels for imregtform, 3 (default)
%       - Interp_Rigid: the interp algorithm when appling transform, ["nearest",
%                   "linear","cubic"(default)], where "bilinear" and "bicubic" is
%                   also valid but as same as "linear" and "cubic"
%       - Interp_NonRigid: the interp algorithm when appling transform, ["nearest",
%                   "linear","cubic"(default)], where "bilinear" and "bicubic" is
%                   also valid but as same as "linear" and "cubic"
%       - LRTDS: the local registration transformation downsampling level,
%                integer number for how many pixel cause a sampling point,
%                1(no downsample) as default
%       - BigFile: the flag for using big file support, please keep it on
%                  if the file size on disk is greater than 10% memory,
%                  value could be "on"/"off", "off"(default)
%   output:
%   - signal: the signal channel has been aligned, 4-D uint16
%   - tform: 1st-col: the global transform object, affine2d or affine3d
%            2nd-col: the local displacement field, double
%   - marker: the marker channel has been aligned, 4-D uint16
%
%   see also: imregdemons, imregconfig, imregtform, imref3d, imregcorr,
%             pca, medfilt3, imboxfilt3, imhistmatchn, affine3d

% REGISTER3D_TC_AUTO:
% Version: 1.0.0
%   *** Basic registration functions
%
% Version: 1.0.1
%   *** the big file (more than physical memory) register is valid
%   *** multi-channel register is supported
%
% Version: 1.0.2
%   *** solve the pyramid levels error when sample small size occured
%   *** filtered movie for more robust when register
%
% Version: 1.0.3
%   *** fix the parallel worker over malloc bug
%
% Version: 1.0.4
%   *** use normcorr2 fix z shift between two volumes
%
% Version: 1.0.5
%   *** use local z direction template matching to speed up estimate
%   *** fix z fixed estimate bug for speed up
%
% Version: 1.0.6
%   *** decrease the memory alloc for single file, remove z estimate helper
%       function for processing speed up
%   *** using parameters checker to make sure inputs is safe and valid
%
% Version: 1.0.7
%   *** fix channel rematch bug
%   *** fix RefVol check bug
%   *** update loadfile usage
%
% Version: 1.0.8
%   *** fix memory allocate bug when cpu register is using
%   *** big file (size over memory) supported
%   *** more register parameters added
%
% Version: 1.0.9
%   *** optimize the registration solution
%   *** auto detect the color channels order, remove the chlOrder selection
%
% Version: 1.1.0
%   *** fix the rigid/similarity moving estimation bug
%   *** append the registration frames series
%
% Version: 1.2.0
%   *** global/local registration parameters full support
%   *** update local registration pre-process
%   *** using LRTDS for decreasing memory allocate when local registration
%
% Version: 1.2.1
%   *** Fix adjustMovie less frames bug
%   *** growth factor is supported for multimodal
%
% Version: 1.2.2
%   *** Fix the big file loading slowly when cut it off
%   *** Fix the pesudo multi-loading bug at unix platform
%   *** Optimize the temporary file loading
%   *** Update the partial registration function
%   *** Fix single frame registration requires parpool bug
%   *** Fix single frame registration padding fail bug
%
% Version: 1.2.3
%   *** Update 'bigfile mode' temporary file redirected user folder 
%   *** Fix Similarity registration error bug on multi-cpu mode
%
% Version: 1.2.4
%   *** Update option: "DS" for downsampling
%   *** Remove border values
% Version: 1.2.5
%   *** Remove and clean up file pointer to a function
%   *** fix file pointer bug

% Copyright (c) 2022-2023, Weihan Li

VALID_CHANNEL_LABEL = ["r", "g", "b"];
VALID_INTERP_METHOD = ["nearest", "linear", "cubic", "bilinear", "bicubic"];
VALID_REGMODE_METHOD = ["global-only", "local-only","mix"];
VALID_MODAL = ["multimodal", "monomodal"];
VALID_TFORM_MODE = ["translation", "rigid", "affine", "similarity"];
VALID_DS_OPTIONS = ["auto","1X1","2X2","3X3"];
VALID_COMP_ACC = [512, 1024, 2048, 4096];

p = inputParser;
valid_file = @(x) isempty(x)||((isstring(x) || ischar(x)) && exist(x,"file"));
valid_movinfo = @(x) isstruct(x) && all(ismember(["mov","opts"],string(fieldnames(x))));
valid_refvol = @(x)isstruct(x) && all(ismember(string(fieldnames(x)),["G","L"]));
valid_regframes = @(x) validateattributes(x,{'numeric'},{'row','positive','integer','increasing'});
valid_regmode = @(x)(ismember(x,VALID_REGMODE_METHOD));
valid_modal = @(x)(ismember(x,VALID_MODAL));
valid_tform = @(x)(ismember(x,VALID_TFORM_MODE));
valid_compacc = @(x)(ismember(x,VALID_COMP_ACC));
valid_initstep = @(x) strcmp(x,"auto") || (isscalar(x) && isreal(x) && x>0);
valid_autocontrast = @(x)validateattributes(x,{'logical'},{'scalar'});
valid_maxitern = @(x)((isvector(x) && isPositiveIntegerValuedNumeric(x))...
    || strcmp(x,"auto"));
valid_minstep = @(x)(isscalar(x) && isreal(x) && x>0);
valid_itercoeff = @(x)validateattributes(x,{'numeric'},{'scalar','real','>',0});
valid_afs = @(x)validateattributes(x,{'numeric'},{'scalar','real','>=',0.5,'<=',3});
valid_coregc = @(x)validateattributes(x,{'numeric'},{'scalar','real','>',0,'<=',100});
valid_ds = @(x)(ismember(x,VALID_DS_OPTIONS));
valid_chla = @(x)isscalar(x) && (isstring(x)||ischar(x)) ...
    && ismember(lower(x),VALID_CHANNEL_LABEL);
valid_chls = @(x)isscalar(x) && (isstring(x)||ischar(x)) ...
    && ismember(lower(x),VALID_CHANNEL_LABEL);
valid_vpl = @(x)validateattributes(x,{'numeric'},...
    {'scalar','integer','>',0,'<',inf});
valid_interp = @(x)(isstring(x)||ischar(x)) && ismember(lower(x),VALID_INTERP_METHOD);
valid_lrtds = @(x)validateattributes(x,{'numeric'},{'scalar','integer'});
valid_bigfile = @(x)(ismember(x,["on","off"]));

%======================== DEFAULT PARAMETER SETTING =======================
default_filename = [];               % the default filename
default_movinfo = struct('mov',[],'opts',[]); % the default movinfo
default_refvol = struct("G",   ["mean", "1"],...
    "L",    ["mean", "1"]);          % default refvol
default_regframes = intmax("uint16");% the default registration frames
default_regmode = "global-only";     % global-only as default, for fast but may not be the most accuracy
default_modal = "monomodal";         % the aligned modal, monomodal as usual
default_tform = "translation";       % translation transformation
default_compacc = 1024;              % the default compensate accuracy
default_initstep = "auto";           % the initial step when optimize registration
default_autocontrast = true;         % the auto contrast flag
default_minstep = 1e-5;              % the minimum step when optimize registration
default_itercoeff = 0.5;             % the iteration coefficient
default_afs = 1.0;                   % the afs
default_max_iter_num = "auto";       % the maximum iteration number when optimize registration
default_co_reg_c = 50;               % the percent of image information when denoised
default_ds = "auto";                 % the default downsampling option
default_channel_align = "r";         % the marker channel (no flashing,weak blanching) need to be aligned
default_channel_signal = "g";        % apply transform on signal channel (flashing,blanching)
default_vpl = 3;                     % the pyramid level number when register
default_interp = "cubic";            % the interp method when applying image operation
default_lrtds = 1;                   % the downsampling coefficient
default_bigfile = "off";             % the big file flag
%==========================================================================

addOptional(p,'movinfo',default_movinfo,valid_movinfo);
addOptional(p,'filename',default_filename,valid_file);
addParameter(p,'RefVol',default_refvol,valid_refvol);
addParameter(p,'RegFrames',default_regframes,valid_regframes);
addParameter(p,'RegMode',default_regmode,valid_regmode);
addParameter(p,'Modal',default_modal,valid_modal);
addParameter(p,'Tform',default_tform,valid_tform);
addParameter(p,'CompAcc',default_compacc,valid_compacc);
addParameter(p,'InitStep',default_initstep,valid_initstep);
addParameter(p,'AutoContrast',default_autocontrast,valid_autocontrast);
addParameter(p,'MinStep',default_minstep,valid_minstep);
addParameter(p,'IterCoeff',default_itercoeff,valid_itercoeff);
addParameter(p,'AFS',default_afs,valid_afs);
addParameter(p,'MaxIterN_Rigid',default_max_iter_num,valid_maxitern);
addParameter(p,'MaxIterN_NonRigid',default_max_iter_num,valid_maxitern);
addParameter(p,'CoRegC',default_co_reg_c,valid_coregc);
addParameter(p,'DS',default_ds,valid_ds);
addParameter(p,'Chla',default_channel_align,valid_chla);
addParameter(p,'Chls',default_channel_signal,valid_chls);
addParameter(p,'VPL_Rigid',default_vpl,valid_vpl);
addParameter(p,'VPL_NonRigid',default_vpl,valid_vpl);
addParameter(p,'Interp_Rigid',default_interp,valid_interp);
addParameter(p,'Interp_NonRigid',default_interp,valid_interp);
addParameter(p,'LRTDS',default_lrtds,valid_lrtds);
addParameter(p,'BigFile',default_bigfile,valid_bigfile);
parse(p,varargin{:});

denoise_dim = [];
% because the different memory management between
% unix(struct with deep copy) and windows(struct with shallow copy)
parser_results = p.Results;

switch parser_results.Modal
    case 'multimodal'
        if isequal(parser_results.InitStep,"auto")
            initstep = 6.25e-3;
        else
            initstep = parser_results.InitStep;
        end
        if isequal(parser_results.IterCoeff, default_itercoeff)
            itercoeff = 1.05;
        else
            itercoeff = parser_results.IterCoeff;
        end
        if isequal(parser_results.MinStep, default_minstep)
            minstep = 1.5e-6;
        else
            minstep = parser_results.MinStep;
        end
    case 'monomodal'
        if isequal(parser_results.InitStep,"auto")
            initstep = 6.25e-2;
        else
            initstep = parser_results.InitStep;
        end
        assert(parser_results.IterCoeff>0 && parser_results.IterCoeff<1,...
            "register3D_auto:invalidInterationCoefficient",...
            "Invalid iteration coefficient for monomodal.");
        itercoeff = parser_results.IterCoeff;
        minstep = parser_results.MinStep;
end

if isequal(parser_results.MaxIterN_Rigid,"auto")
    switch parser_results.Modal
        case 'multimodal'
            max_iter_num_rigid = round(100*6.25e-3/initstep);
        case 'monomodal'
            max_iter_num_rigid = round(100*6.25e-2/initstep);
    end
else
    max_iter_num_rigid = parser_results.MaxIterN_Rigid;
end

if isequal(parser_results.MaxIterN_NonRigid, "auto")
    max_iter_num_nonrigid = 100;
else
    max_iter_num_nonrigid = parser_results.MaxIterN_NonRigid;
end

% load the file
if isempty(parser_results.movinfo.mov)
    [opts,mov,~,~] = loadfile(parser_results.filename);
else
    mov = parser_results.movinfo.mov;
    opts = parser_results.movinfo.opts;
end

channel_order = opts.cOrder;

% generate global value: fixed volume and data scale
[MA_MIN, MA_MAX] = getMinMaxIn(mov, parser_results.Chla, channel_order);
fv = GenFixVolProfile(mov, parser_results, channel_order, parser_results.DS);

% pre-processing data that was RegFrames mentioned
if isequal(parser_results.RegFrames,intmax("uint16"))
    reg_frames = 1:opts.frames;
else
    reg_frames = intersect(parser_results.RegFrames, 1:opts.frames);
end
mov(:,:,:,:,setdiff(1:opts.frames,reg_frames)) = []; % remove the kept data
opts.frames = numel(reg_frames); % adjust the frames
opts.images = opts.frames*opts.slices*opts.channels; % adjust the total images

CONTRAST_CONSTANT = parser_results.CompAcc;

reg_param = struct( 'regmode',      parser_results.RegMode,...
    'fixvol',        fv,...
    'modal',        parser_results.Modal,...
    'tform',        parser_results.Tform,...
    'initStep',     initstep,...
    'autoContrast', parser_results.AutoContrast,...
    'minStep',      minstep,...
    'iterCoeff',    itercoeff,...
    'afs',          parser_results.AFS,...
    'maxIterNumRigid',   max_iter_num_rigid,...
    'maxIterNumNonRigid',max_iter_num_nonrigid,...
    'coRegC',       parser_results.CoRegC,...
    'DS',           parser_results.DS,...
    'chlA',         parser_results.Chla,...
    'chlS',         parser_results.Chls,...
    'chlMode',      channel_order,...
    'vPLRigid',     parser_results.VPL_Rigid,...
    'vPLNonRigid',  parser_results.VPL_NonRigid,...
    'interpRigid',  parser_results.Interp_Rigid,...
    'interpNonRigid',parser_results.Interp_NonRigid,...
    'LRTDS',        parser_results.LRTDS,...
    'bigfile',      parser_results.BigFile);

% clear vars for debug easy
clearvars -except opts mov reg_param CONTRAST_CONSTANT MA_MIN MA_MAX;

if ~isempty(mov)
    switch reg_param.bigfile
        case "on"
            chl = [find(reg_param.chlA==reg_param.chlMode), ...
                   find(reg_param.chlS==reg_param.chlMode)];
            [ma_ptr, ms_ptr] = GenFilePointer(mov, chl);
        case "off"
            % use the memory variable directly
            ma_ptr = squeeze(uint16(mov(:,:,...
                reg_param.chlA==reg_param.chlMode,:,:)));
            ms_ptr = squeeze(uint16(mov(:,:,...
                reg_param.chlS==reg_param.chlMode,:,:)));
        otherwise
    end
    clearvars mov;  % free memory for more computation

    [tform,signal,marker] = reg(ma_ptr,ms_ptr,opts,reg_param);
else
    signal = [];
    marker = [];
    tform = [];
    return;
end

    function [tform,ms_ptr,ma_ptr] = reg(ma_ptr,ms_ptr,opts,reg_param)
        % input:
        %   - ma_ptr: the movie marker channel (or disk file pointer)
        %   - ms_ptr: the movie signal channel (or disk file pointer)
        %   - opts:   the options of movie file
        %   - reg_param: the registration paramaters
        % output:
        %   - tform:  the geometric transformation
        %   - ms_ptr: the aligned movie signal channel (or disk file pointer)
        %   - ma_ptr: the aligned movie marker channel (or disk file pointer)

        % extract the fixvol
        fixvol_global_a_ds = reg_param.fixvol.fixvol_global_a_ds;
        ds_scale = reg_param.fixvol.ds_scale;
        fixvol_local_a = reg_param.fixvol.fixvol_local_a;
        borderval_global_a = reg_param.fixvol.borderval_global_a;
        borderval_global_s = reg_param.fixvol.borderval_global_s;

        tform = cell(opts.frames,2);

        [rm,cn,~] = getRunningMode();

        % !!! HERE IS THE FIRST DENOISE USING !!!
        [fixvolEst_ds,denoise_dim] = denoise3D(fixvol_global_a_ds,'CRC',reg_param.coRegC);

        Reg3D = imref3d(size(fixvol_local_a),opts.xRes,opts.yRes,opts.zRes);
        Reg3D_ds  = imref3d(size(fixvol_global_a_ds),opts.xRes/ds_scale,...
            opts.yRes/ds_scale, opts.zRes/1);

        % turn off the imregcorr 'weakpeak' warning
        warning('off','images:imregcorr:weakPeakCorrelation');

        % select in difference register mode
        % NOTE: YOU CAN HARDLY GO INTO 'CPU' MODE, JUST FOR DEBUGING
        % REGISTRATION WITH NON-RIGID WILL MOFDIFY THE ALIGN CHANNEL JUST
        % IN THE 'CPU' MODE
        if isequal(rm,'cpu') || opts.frames == 1
            switch reg_param.bigfile
                case "on"
                    warning("Are you crazy? " + ...
                        "Why do you process such a big file with single core cpu? " + ...
                        "Please, DO NOT PUA your computer any more.");
                    return;
                case "off"
                    % set a waitbar (you can cancel process immediately)
                    bar = waitbar(0,"start to register volume...","Name","Relax bar",...
                        "CreateCancelBtn",'setappdata(gcbf,''canceling'',true)');
                    setappdata(bar,'canceling',false);

                    % Set up the Initial Registration
                    [optimizer, metric] = imregconfig(reg_param.modal);
                    if isequal(reg_param.modal,"multimodal")
                        optimizer.InitialRadius = reg_param.initStep;
                        optimizer.Epsilon = reg_param.minStep;
                        optimizer.GrowthFactor = reg_param.iterCoeff;
                    else
                        optimizer.MaximumStepLength = reg_param.initStep;
                        optimizer.MinimumStepLength = reg_param.minStep;
                        optimizer.RelaxationFactor = reg_param.iterCoeff;
                    end
                    optimizer.MaximumIterations = reg_param.maxIterNumRigid;

                    % genearte the downsampling data
                    ma_ptr_ds = DownSampling(ma_ptr, ds_scale);

                    % align each volume with the fixed volume
                    for k = 1:opts.frames
                        if getappdata(bar,'canceling')
                            break;
                        end

                        switch reg_param.regmode
                            case 'global-only'
                                %%%%%%%%% global align %%%%%%%%
                                % coarse registration transformation -> (x,y) translation
                                % [WARNING] tform3D is affine3d object
                                if reg_param.tform == "translation"
                                    pseudo_tform3D_ds = coreg3D(ma_ptr_ds(:,:,:,k),...
                                        Reg3D_ds, fixvolEst_ds, borderval_global_a,...
                                        denoise_dim, reg_param, opts);

                                    tform{k,1} = imregtform(ma_ptr_ds(:,:,:,k), Reg3D_ds, fixvol_global_a_ds,...
                                        Reg3D_ds, reg_param.tform, optimizer, metric, "PyramidLevels", reg_param.vPLRigid,...
                                        "InitialTransformation",pseudo_tform3D_ds);
                                else
                                    tform{k,1} = imregtform(ma_ptr_ds(:,:,:,k), Reg3D_ds, fixvol_global_a_ds,...
                                        Reg3D_ds, reg_param.tform, optimizer, metric, "PyramidLevels", reg_param.vPLRigid);
                                end

                                ma_ptr(:,:,:,k) = imwarp(ma_ptr(:,:,:,k),Reg3D, ...
                                    tform{k,1},reg_param.interpRigid,'OutputView',Reg3D,...
                                    'SmoothEdge',true,'FillValues',borderval_global_a);
                            case 'local-only'
                                % pre-process
                                ma_ptr(:,:,:,k) = imrescalei(ma_ptr(:,:,:,k),...
                                    MA_MIN, MA_MAX, true);

                                ma_ptr(:,:,:,k) = imhistmatchn(ma_ptr(:,:,:,k),...
                                    fixvol_local_a, CONTRAST_CONSTANT);

                                % local align
                                [tform{k,2},~] = imregdemons(ma_ptr(:,:,:,k),...
                                    fixvol_local_a, reg_param.maxIterNumNonRigid,...
                                    "AccumulatedFieldSmoothing",reg_param.afs,...
                                    "PyramidLevels",reg_param.vPLNonRigid);

                                % post-process, recover distribution
                                ma_ptr(:,:,:,k) = imrescalei(ma_ptr(:,:,:,k),...
                                    MA_MIN, MA_MAX, false);

                                ma_ptr(:,:,:,k) = imwarp(ma_ptr(:,:,:,k),tform{k,2},...
                                    reg_param.interpNonRigid,'SmoothEdge',true);
                            case 'mix'
                                % [WARNING] tform3D is affine3d object
                                if reg_param.tform == "translation"
                                    % pseudo tform3D for translation and affine
                                    pseudo_tform3D_ds = coreg3D(ma_ptr_ds(:,:,:,k),...
                                        Reg3D_ds, fixvolEst_ds, borderval_global_a,...
                                        denoise_dim, reg_param, opts);

                                    tform{k,1} = imregtform(ma_ptr_ds(:,:,:,k),Reg3D_ds,fixvol_global_a_ds,...
                                        Reg3D_ds,reg_param.tform,optimizer,metric,"PyramidLevels",reg_param.vPLRigid,...
                                        "InitialTransformation",pseudo_tform3D_ds);
                                else
                                    tform{k,1} = imregtform(ma_ptr_ds(:,:,:,k),Reg3D_ds,fixvol_global_a_ds,...
                                        Reg3D_ds,reg_param.tform,optimizer,metric,"PyramidLevels",reg_param.vPLRigid);
                                end
                                ma_ptr(:,:,:,k) = imwarp(ma_ptr(:,:,:,k),Reg3D, tform{k,1},...
                                    reg_param.interpRigid,'OutputView',Reg3D,'SmoothEdge',true,...
                                    'FillValues',borderval_global_a);

                                % pre-process
                                ma_ptr(:,:,:,k) = imrescalei(ma_ptr(:,:,:,k),...
                                    MA_MIN, MA_MAX, true);

                                ma_ptr(:,:,:,k) = imhistmatchn(ma_ptr(:,:,:,k),...
                                    fixvol_local_a, CONTRAST_CONSTANT);

                                [tform{k,2},~] = imregdemons(ma_ptr(:,:,:,k),...
                                    fixvol_local_a, reg_param.maxIterNumNonRigid,...
                                    "AccumulatedFieldSmoothing",reg_param.afs,...
                                    "PyramidLevels",reg_param.vPLNonRigid);

                                % post-process
                                ma_ptr(:,:,:,k) = imrescalei(ma_ptr(:,:,:,k),...
                                    MA_MIN, MA_MAX, false);

                                ma_ptr(:,:,:,k) = imwarp(ma_ptr(:,:,:,k),tform{k,2},...
                                    reg_param.interpNonRigid,'SmoothEdge',true);
                            otherwise
                                error("Unsupported registration mode.");
                        end

                        % apply transform on mov_signal
                        if ~isempty(tform{k,1})
                            ms_ptr(:,:,:,k) = imwarp(ms_ptr(:,:,:,k),Reg3D,...
                                tform{k,1},reg_param.interpRigid,...
                                'OutputView',Reg3D,'SmoothEdge',true,...
                                'FillValues', borderval_global_s);
                        end
                        if ~isempty(tform{k,2})
                            ms_ptr(:,:,:,k) = imwarp(ms_ptr(:,:,:,k),...
                                tform{k,2},reg_param.interpNonRigid,'SmoothEdge',true);
                        end

                        % transform from double to single and downsampling
                        % for decreasing memory allocation
                        tform{k,2} = DownsamplingDisplacement(tform{k,2}, reg_param.LRTDS);

                        waitbar(k/opts.frames,bar,"register processing "...
                            +num2str(k/opts.frames*100,3)+"% ...");
                    end

                    if k == opts.frames,msg = "register succeed.";else, ...
                            msg = "register canceled.";end
                    waitbar(k/opts.frames,bar,msg);
                    pause(1);
                    delete(bar);
                otherwise
            end
        else
            switch reg_param.regmode
                case 'global-only'
                    bar = parwaitbar(opts.frames,'Waitbar',true);
                    [tform,ms_ptr,ma_ptr] = align_global_cpu(ma_ptr,...
                        ms_ptr,tform);
                case 'local-only'
                    bar = parwaitbar(opts.frames,'Waitbar',true);
                    if isequal(rm,'gpu') || isequal(rm,'multi-gpu')
                        [tform,ms_ptr,ma_ptr] = align_local_gpu(ma_ptr,...
                            ms_ptr,tform);
                    else
                        [tform,ms_ptr,ma_ptr] = align_local_cpu(ma_ptr,...
                            ms_ptr,tform);
                    end
                case 'mix'
                    bar = parwaitbar(2*opts.frames,'Waitbar',true);

                    [tform,ms_ptr,ma_ptr] = align_global_cpu(ma_ptr,...
                        ms_ptr,tform);

                    if isequal(rm,'gpu') || isequal(rm,'multi-gpu')
                        [tform,ms_ptr,ma_ptr] = align_local_gpu(ma_ptr,...
                            ms_ptr,tform);
                    else
                        [tform,ms_ptr,ma_ptr] = align_local_cpu(ma_ptr,...
                            ms_ptr,tform);
                    end
                otherwise
                    error("unsupported register mode");
            end
            bar.Destroy;

            % turn on the imregcorr 'weakpeak' warning
            warning('on','images:imregcorr:weakPeakCorrelation');

            % clear Functions: clears all Parallel Computing Toolbox objects from the current MATLAB session
            clear functions %#ok<CLFUNC> 
        end

        function [t,ms_new,ma_new] = align_global_cpu(ma,ms,t)
            parobj = OpenParpool(min(min(cn,getCpuBlockNumber(reg_param.bigfile)),opts.frames));

            mpiprofile on

            if reg_param.bigfile == "on"
                % using for loop load block data
                mv_size = size(ma,'mov_aligned');
                % generate the registration results as preprocess
                folder_local = fileparts(ma.Properties.Source);
                ms_new = matfile(fullfile(folder_local, 'mov_signal_global.mat'), ...
                    "Writable",true);
                ms_new.mov_signal(mv_size(1),mv_size(2),mv_size(3),mv_size(4))...
                    = uint16(0);        % system malloc standby
                ma_new = matfile(fullfile(folder_local, 'mov_aligned_global.mat'), ...
                    "Writable",true);
                ma_new.mov_aligned(mv_size(1),mv_size(2),mv_size(3),mv_size(4))...
                    = uint16(0);
                NWorker = parobj.NumWorkers;
                NB = ceil(mv_size(end)/NWorker);
                for ii = 1:NB
                    bs = min(NWorker,mv_size(end)-(ii-1)*NWorker);
                    ma_block = ma.mov_aligned(:,:,:,...
                        (ii-1)*NWorker+1:(ii-1)*NWorker+bs);
                    ma_block_ds = DownSampling(ma_block, ds_scale);
                    ms_block = ms.mov_signal(:,:,:,...
                        (ii-1)*NWorker+1:(ii-1)*NWorker+bs);
                    t_block = cell(bs,1);
                    parfor m = 1:bs
                        % Set up the Initial Registration
                        [optimizer_par,metric_par] = imregconfig(reg_param.modal); %#ok<PFBNS>
                        if isequal(reg_param.modal,"multimodal")
                            optimizer_par.InitialRadius = reg_param.initStep;
                            optimizer_par.Epsilon = reg_param.minStep;
                            optimizer_par.GrowthFactor = reg_param.iterCoeff;
                        elseif isequal(reg_param.modal,"monomodal")
                            optimizer_par.MaximumStepLength = reg_param.initStep;
                            optimizer_par.MinimumStepLength = reg_param.minStep;
                            optimizer_par.RelaxationFactor = reg_param.iterCoeff;
                        end
                        optimizer_par.MaximumIterations = reg_param.maxIterNumRigid;

                        if reg_param.tform == "translation"
                            % TMPDATA MUST BE CLAIM EXPLICITLY IN PARFOR!!!
                            % coarse registration
                            pstf3D_ds = coreg3D(ma_block_ds(:,:,:,m),...
                                Reg3D_ds, fixvolEst_ds, borderval_global_a, ...
                                denoise_dim, reg_param, opts);

                            % fine registration
                            tformk1 = imregtform(ma_block_ds(:,:,:,m),Reg3D_ds,fixvol_global_a_ds,...
                                Reg3D_ds,reg_param.tform,optimizer_par,metric_par,...
                                "PyramidLevels",reg_param.vPLRigid,...
                                "InitialTransformation",pstf3D_ds);
                        else
                            % fine registration
                            tformk1 = imregtform(ma_block_ds(:,:,:,m),Reg3D_ds,fixvol_global_a_ds,...
                                Reg3D_ds,reg_param.tform,optimizer_par,metric_par,...
                                "PyramidLevels",reg_param.vPLRigid);
                        end
                        
                        % warp the geometric transformation
                        ms_block(:,:,:,m) = imwarp(ms_block(:,:,:,m),...
                            Reg3D,tformk1,reg_param.interpRigid,...
                            'OutputView',Reg3D,'SmoothEdge',true,...
                            'FillValues', borderval_global_s);
                        ma_block(:,:,:,m) = imwarp(ma_block(:,:,:,m),...
                            Reg3D,tformk1,reg_param.interpRigid,...
                            'OutputView',Reg3D,'SmoothEdge',true,...
                            'FillValues',borderval_global_a);
                        t_block(m,1) = {tformk1};
                        bar.Send; %#ok<PFBNS>
                    end
                    ms_new.mov_signal(:,:,:,(ii-1)*NWorker+1:(ii-1)*NWorker+bs)...
                        = ms_block;
                    ma_new.mov_aligned(:,:,:,(ii-1)*NWorker+1:(ii-1)*NWorker+bs)...
                        = ma_block;
                    t((ii-1)*NWorker+1:(ii-1)*NWorker+bs) = t_block;
                end
            else
                parfor m = 1:opts.frames
                    % Set up the Initial Registration
                    [optimizer_par,metric_par] = imregconfig(reg_param.modal); %#ok<PFBNS>
                    if isequal(reg_param.modal,"multimodal")
                        optimizer_par.InitialRadius = reg_param.initStep;
                        optimizer_par.Epsilon = reg_param.minStep;
                        optimizer_par.GrowthFactor = reg_param.iterCoeff;
                    elseif isequal(reg_param.modal,"monomodal")
                        optimizer_par.MaximumStepLength = reg_param.initStep;
                        optimizer_par.MinimumStepLength = reg_param.minStep;
                        optimizer_par.RelaxationFactor = reg_param.iterCoeff;
                    end
                    optimizer_par.MaximumIterations = reg_param.maxIterNumRigid;

                    % genearte the down sampling aligned channel data
                    ma_m_ds = DownSampling(ma(:,:,:,m), ds_scale);

                    if reg_param.tform == "translation"
                        % TMPDATA MUST BE CLAIM EXPLICITLY IN PARFOR!!!
                        % coarse registration
                        pstf3D_ds = coreg3D(ma_m_ds,...
                            Reg3D_ds, fixvolEst_ds, borderval_global_a, ...
                            denoise_dim, reg_param, opts);

                        % fine registration
                        tformk1 = imregtform(ma_m_ds,Reg3D_ds,fixvol_global_a_ds,...
                            Reg3D_ds,reg_param.tform,optimizer_par,metric_par,...
                            "PyramidLevels",reg_param.vPLRigid,...
                            "InitialTransformation",pstf3D_ds);
                    else
                        % fine registration
                        tformk1 = imregtform(ma_m_ds,Reg3D_ds,fixvol_global_a_ds,...
                            Reg3D_ds,reg_param.tform,optimizer_par,metric_par,...
                            "PyramidLevels",reg_param.vPLRigid);
                    end
                    
                    ms(:,:,:,m) = imwarp(ms(:,:,:,m),Reg3D,tformk1,...
                        reg_param.interpRigid,'OutputView',Reg3D,...
                        'SmoothEdge',true,'FillValues',borderval_global_s);
                    ma(:,:,:,m) = imwarp(ma(:,:,:,m),Reg3D,tformk1,...
                        reg_param.interpRigid,'OutputView',Reg3D,...
                        'SmoothEdge',true,'FillValues',borderval_global_a);
                    t(m,1) = {tformk1};
                    bar.Send; %#ok<PFBNS>
                end
                % please do not deep copy, linux! thank you so much!!!
                ms_new = ms;
                ma_new = ma;
            end

            mpiprofile viewer

            CloseParpool(parobj);
        end

        function [t,ms_new,ma_new] = align_local_gpu(ma,ms,t)
            parobj = OpenParpool(min(getGpuBlockNumber(),opts.frames));

            mpiprofile on

            if reg_param.bigfile == "on"
                % using for loop load block data
                mv_size = size(ma,'mov_aligned');
                % generate the registration results as preprocess
                folder_local = fileparts(ma.Properties.Source);
                ms_new = matfile(fullfile(folder_local, 'mov_signal_local.mat'), ...
                    "Writable",true);
                ms_new.mov_signal(mv_size(1),mv_size(2),mv_size(3),mv_size(4))...
                    = uint16(0);        % system malloc standby
                ma_new = matfile(fullfile(folder_local, 'mov_aligned_local.mat'), ...
                    "Writable",true);
                ma_new.mov_aligned(mv_size(1),mv_size(2),mv_size(3),mv_size(4))...
                    = uint16(0);
                NWorker = parobj.NumWorkers;
                NB = ceil(mv_size(end)/NWorker);
                for ii = 1:NB
                    bs = min(NWorker,mv_size(end)-(ii-1)*NWorker);

                    ma_block = ma.mov_aligned(:,:,:,...
                        (ii-1)*NWorker+1:(ii-1)*NWorker+bs);

                    % pre-process
                    ma_block = imrescalei(ma_block, MA_MIN, MA_MAX, true);
                    if reg_param.autoContrast == true
                        for m = 1:bs
                            ma_block(:,:,:,m) = imhistmatchn(ma_block(:,:,:,m),...
                                fixvol_local_a, CONTRAST_CONSTANT);
                        end
                    end

                    ms_block = ms.mov_signal(:,:,:,...
                        (ii-1)*NWorker+1:(ii-1)*NWorker+bs);

                    t_block = cell(bs,1);

                    parfor m = 1:bs
                        fixvol_ga = gpuArray(fixvol_local_a);
                        movie_ga = gpuArray(ma_block(:,:,:,m));

                        % GPU running
                        [tformk2,~] = imregdemons(movie_ga,fixvol_ga, ...
                            reg_param.maxIterNumNonRigid,...
                            "AccumulatedFieldSmoothing",reg_param.afs,...
                            "PyramidLevels",reg_param.vPLNonRigid);  %#ok<PFBNS>
                        t_block(m,1) = {gather(tformk2)};
                        % CPU running
                        ma_block(:,:,:,m) = imwarp(ma_block(:,:,:,m),t_block{m,1},...
                            reg_param.interpNonRigid,'SmoothEdge',true);
                        ms_block(:,:,:,m) = imwarp(ms_block(:,:,:,m),t_block{m,1},...
                            reg_param.interpNonRigid,'SmoothEdge',true);

                        bar.Send; %#ok<PFBNS>
                    end

                    % post-process 1
                    ma_block = imrescalei(ma_block, MA_MIN, MA_MAX, false);

                    ms_new.mov_signal(:,:,:,(ii-1)*NWorker+1:(ii-1)*NWorker+bs)...
                        = ms_block;
                    ma_new.mov_aligned(:,:,:,(ii-1)*NWorker+1:(ii-1)*NWorker+bs)...
                        = ma_block;

                    % post-process 2
                    for m = 1:bs
                        t{(ii-1)*NWorker + m, 2} = DownsamplingDisplacement(t_block{m},...
                            reg_param.LRTDS);
                    end
                end
            else
                %pre-process
                ma = imrescalei(ma, MA_MIN, MA_MAX, true);

                parfor m = 1:opts.frames
                    % balance CPU-GPU interection
                    % 1. GPU register, CPU imwarp with linear, fast but low accuracy
                    % 2.[default] GPU register, GPU imwarp with cubic, slow but high accuracy
                    % [MAY CAUSING BUG] GPU makes align channel interp method
                    % is 'linear', but CPU makes signal channel interp method is reg_param.interp
                    % send array to GPU
                    fixvol_ga = gpuArray(fixvol_local_a);

                    if reg_param.autoContrast == true %#ok<PFBNS>
                        ma(:,:,:,m) = imhistmatchn(ma(:,:,:,m),...
                            fixvol_local_a, CONTRAST_CONSTANT);
                    end
                    movie_ga = gpuArray(ma(:,:,:,m));

                    % GPU running
                    [tformk2,~] = imregdemons(movie_ga,fixvol_ga, ...
                        reg_param.maxIterNumNonRigid,...
                        "AccumulatedFieldSmoothing",reg_param.afs,...
                        "PyramidLevels",reg_param.vPLNonRigid);
                    t(m,2) = {gather(tformk2)};
                    % CPU running
                    ma(:,:,:,m) = imwarp(ma(:,:,:,m),t{m,2},reg_param.interpNonRigid,...
                        'SmoothEdge',true);
                    ms(:,:,:,m) = imwarp(ms(:,:,:,m),t{m,2},reg_param.interpNonRigid,...
                        'SmoothEdge',true);

                    bar.Send; %#ok<PFBNS>
                end

                % post-process 1
                ma = imrescalei(ma, MA_MIN, MA_MAX, false);

                % post-process 2
                for m = 1:opts.frames
                    t{m, 2} = DownsamplingDisplacement(t{m, 2}, reg_param.LRTDS);
                end

                ms_new = ms;
                ma_new = ma;
            end

            mpiprofile viewer

            CloseParpool(parobj);
        end

        function [t,ms_new,ma_new] = align_local_cpu(ma,ms,t)
            parobj = OpenParpool(min(min(cn,getCpuBlockNumber(reg_param.bigfile)),opts.frames));

            mpiprofile on

            if reg_param.bigfile == "on"
                % using for loop load block data
                mv_size = size(ma,'mov_aligned');
                % generate the registration results as preprocess
                folder_local = fileparts(ma.Properties.Source);
                ms_new = matfile(fullfile(folder_local, 'mov_signal_local.mat'), ...
                    "Writable",true);
                ms_new.mov_signal(mv_size(1),mv_size(2),mv_size(3),mv_size(4))...
                    = uint16(0);        % system malloc standby
                ma_new = matfile(fullfile(folder_local, 'mov_aligned_local.mat'), ...
                    "Writable",true);
                ma_new.mov_aligned(mv_size(1),mv_size(2),mv_size(3),mv_size(4))...
                    = uint16(0);
                NWorker = parobj.NumWorkers;
                NB = ceil(mv_size(end)/NWorker);
                for ii = 1:NB
                    bs = min(NWorker,mv_size(end)-(ii-1)*NWorker);
                    ma_block = ma.mov_aligned(:,:,:,...
                        (ii-1)*NWorker+1:(ii-1)*NWorker+bs);

                    % precess the ma outer loop
                    ma_block = imrescalei(ma_block, MA_MIN, MA_MAX, true);
                    if reg_param.autoContrast == true
                        for m = 1:bs
                            ma_block(:,:,:,m) = imhistmatchn(ma_block(:,:,:,m),...
                                fixvol_local_a,CONTRAST_CONSTANT);
                        end
                    end
                    ms_block = ms.mov_signal(:,:,:,...
                        (ii-1)*NWorker+1:(ii-1)*NWorker+bs);
                    t_block = cell(bs,1);
                    parfor m = 1:bs
                        [tformk2,~] = imregdemons(ma_block(:,:,:,m),fixvol_local_a,...
                            reg_param.maxIterNumNonRigid,...
                            'DisplayWaitbar',false,...
                            'AccumulatedFieldSmoothing',reg_param.afs,...
                            "PyramidLevels",reg_param.vPLNonRigid); %#ok<PFBNS>
                        t_block(m,1) = {tformk2};
                        ma_block(:,:,:,m) = imwarp(ma_block(:,:,:,m),...
                            tformk2,reg_param.interpNonRigid,'SmoothEdge',true);
                        ms_block(:,:,:,m) = imwarp(ms_block(:,:,:,m),...
                            tformk2,reg_param.interpNonRigid,'SmoothEdge',true);

                        bar.Send; %#ok<PFBNS>
                    end
                    % post-process 1
                    ma_block = imrescalei(ma_block, MA_MIN, MA_MAX, false);

                    ms_new.mov_signal(:,:,:,(ii-1)*NWorker+1:(ii-1)*NWorker+bs)...
                        = ms_block;
                    ma_new.mov_aligned(:,:,:,(ii-1)*NWorker+1:(ii-1)*NWorker+bs)...
                        = ma_block;

                    % post-process 2
                    for m = 1:bs
                        t{(ii-1)*NWorker + m, 2} = DownsamplingDisplacement(t_block{m},...
                            reg_param.LRTDS);
                    end
                end
            else
                % pre-process
                ma = imrescalei(ma, MA_MIN, MA_MAX, true);

                parfor m = 1:opts.frames
                    if reg_param.autoContrast == true %#ok<PFBNS>
                        ma(:,:,:,m) = imhistmatchn(ma(:,:,:,m),...
                            fixvol_local_a, CONTRAST_CONSTANT);
                    end
                    [tformk2,~] = imregdemons(ma(:,:,:,m),fixvol_local_a,...
                        reg_param.maxIterNumNonRigid,...
                        'DisplayWaitbar',false,...
                        'AccumulatedFieldSmoothing',reg_param.afs,...
                        "PyramidLevels",reg_param.vPLNonRigid);
                    t(m,2) = {tformk2};
                    ma(:,:,:,m) = imwarp(ma(:,:,:,m),tformk2,reg_param.interpNonRigid,...
                        'SmoothEdge',true);
                    ms(:,:,:,m) = imwarp(ms(:,:,:,m),tformk2,reg_param.interpNonRigid,...
                        'SmoothEdge',true);

                    bar.Send; %#ok<PFBNS>
                end

                % post-process 1
                ma = imrescalei(ma, MA_MIN, MA_MAX, false);

                % post-process 2
                for m = 1:opts.frames
                    t{m, 2} = DownsamplingDisplacement(t{m, 2}, reg_param.LRTDS);
                end

                ms_new = ms;
                ma_new = ma;
            end

            mpiprofile viewer

            CloseParpool(parobj);
        end
    end

    function [MA_MIN, MA_MAX] = getMinMaxIn(mov, ch_slt, ch_ord)
        MA_MIN = min(mov(:,:,ch_slt==ch_ord,:,:),[],"all");
        MA_MAX = max(mov(:,:,ch_slt==ch_ord,:,:),[],"all");
    end

    function fv = GenFixVolProfile(mov, ps, cd, ds)
        tmp_ma = squeeze(mov(:,:,ps.Chla==cd,:,:));
        fixvol_global_a = GetFixedVol(tmp_ma, ps.RefVol.G);
        fixvol_local_a = GetFixedVol(tmp_ma, ps.RefVol.L);
        clearvars tmp_ma;

        tmp_ms = squeeze(mov(:,:,ps.Chls==cd,:,:));
        fixvol_global_s = GetFixedVol(tmp_ms, ps.RefVol.G);
        clearvars tmp_ms;
        
        fixvol_local_a = imrescalei(fixvol_local_a, MA_MIN, MA_MAX, true);

        % get the image border fill value =========================
        borderval_global_a = imborderval(fixvol_global_a, [5,2,5,2,5,5]);
        borderval_global_a = mean(borderval_global_a([1,3,5,6]));
        borderval_global_s = imborderval(fixvol_global_s, [5,2,5,2,5,5]);
        borderval_global_s = mean(borderval_global_s([1,3,5,6]));

        fv.fixvol_local_a = fixvol_local_a;
        fv.borderval_global_a = borderval_global_a;
        fv.borderval_global_s = borderval_global_s;

        switch ds
            case "auto"
                [fv.fixvol_global_a_ds, fv.ds_scale] = DownSampling(fixvol_global_a);
            otherwise
                [fv.fixvol_global_a_ds, fv.ds_scale] = DownSampling(...
                    fixvol_global_a, 1/str2double(ds.extract(1)));
        end
    end

% this function for detecting the best running mode for local hardware setting
    function [rm, workers_n_stable, gpu_cores_n] = getRunningMode()
        % detect the hardware and output the best section
        cpu_cores_n=feature('numCores');
        if cpu_cores_n <= 8
            workers_n_stable = cpu_cores_n;
        else
            % TODO:
            % worker_n_stable depends on image file size
            workers_n_stable = round(0.85*cpu_cores_n);
        end
        gpu_cores_n = gpuDeviceCount('available');
        fprintf('-> available cpu core: %d \n-> available gpu: %d\n',...
            cpu_cores_n,gpu_cores_n);
        fprintf('-> maximum stable workers number: %d\n',workers_n_stable);
        if gpu_cores_n > 0
            if gpu_cores_n == 1
                rm = 'gpu';
            else
                rm = 'multi-gpu';
            end
        else
            if cpu_cores_n > 1
                rm = 'multi-cpu';
            else
                % you know, single core for debugging forever
                rm = 'cpu';
            end
        end
        % Next line for debugging
        % rm = 'cpu';
    end

    function gn = getGpuBlockNumber()
        % This function get the avaiable gpu memory size and give a bench
        % size of a processing job, which is the frame could be loaded
        SINGLE_BYTES = 4;
        FOLD_RATIO = 160;   % the linear estimation may be wrong
        SECURATY_RATIO = 0.85;

        % THE MAXSIZE OF GRAPHICS CARD MEMORY CONTAINS MODEL CAN BE
        % CALCULATE BY LINEAR SIMILARITY

        gpu_memory = zeros(1,gpuDeviceCount);
        for k = 1:gpuDeviceCount
            c = gpuDevice(k);
            gpu_memory(k) = c.AvailableMemory;
        end
        minimal_aval_memory = min(gpu_memory);

        mem_per_volume = opts.width*opts.height*opts.slices ...
            *SINGLE_BYTES;

        if ispc()
            mem_per_worker = 800*1024*1024; % bytes
        elseif isunix()
            mem_per_worker = 1100*1024*1024; %bytes
        else
            % what is the fuck?
            error('Your operation system is so coooool.');
        end

        gn = round(gpuDeviceCount*minimal_aval_memory*SECURATY_RATIO/...
            (mem_per_volume*FOLD_RATIO+mem_per_worker));
        % for even gn to balance
        if gn < 1 || mem_per_volume > minimal_aval_memory
            error('No available enough gpu memory.');
        else
            if gn >gpuDeviceCount
                % balance different GPU loading
                gn = gn - mod(gn,gpuDeviceCount);
            end
        end

    end

% TODO: fold ratio relys on whether the data is big file
    function cn = getCpuBlockNumber(useBigFile)
        if ~exist("isBigFile","var") || isempty(useBigFile)
            useBigFile = "off";
        elseif ~isstring(useBigFile)
            error("Invalid input argument.");
        end
        % This function get the avaiable memory size and give a bench
        % size of a processing job, which is the frame could be loaded
        UINT16_BYTES = 2;

        if useBigFile == "on"
            FOLD_RATIO = 2+3;
        elseif useBigFile == "off"
            % 2 for two channel(aligned,signal),2 for argument and temporary
            % variable (pyramid levels algorithm), 3 for tmpdata (coreg3D)
            FOLD_RATIO = 2*2+3;
        end

        % left ?% avaiable memory for system running AND O(1) memory
        % alloc
        SECURATY_RATIO = 0.8;

        % THE MAXSIZE OF MEMORY CONTAINS MODEL CAN BE
        % CALCULATE BY LINEAR SIMILARITY
        mem_per_volume = opts.width*opts.height*opts.slices ...
            *UINT16_BYTES;

        % SELECT PRESENT PLATFORM
        if ispc()
            mem_per_worker = 800*1024*1024; % bytes
            cn = fix(GetAvailableMemory()*SECURATY_RATIO...
                /(mem_per_volume*FOLD_RATIO+mem_per_worker));
        elseif isunix()
            mem_per_worker = 1100*1024*1024; %bytes
            cn = fix(GetAvailableMemory()*SECURATY_RATIO...
                /(mem_per_volume*FOLD_RATIO+mem_per_worker));
        else
            % what is the fuck?
            error('Your operation system is so coooool.');
        end
        if cn < 1
            error('No available enough memory.');
        end
    end

end

% About maximum 3 fold memory of movol used <coreg3D>
function [tformEst,movolEst] = coreg3D(movol,Reg3D,fixvol,bdv,denoise_dim,reg_param,opts)
% this function register movol to fixvol by using projection
% phase corrleation method
fixed_plane = zProjection(fixvol);

Rfixed2D = imref2d(size(fixed_plane),opts.xRes,opts.yRes);

% 1 fold memory copy from <denoise3D>. 2 fold memory alloc when achieving
% peak
movol_denoised = denoise3D(movol,'dd',denoise_dim,'CRC',reg_param.coRegC);
moving_plane = zProjection(movol_denoised);

Rmoving2D = imref2d(size(moving_plane),opts.xRes,opts.yRes);
if isequal(reg_param.tform,"affine")
    tf2D = 'similarity';    % affine, similarity -> similarity
elseif isequal(reg_param.tform,"translation")
    tf2D = reg_param.tform; %
else
    error("register3D:coreg3D Bad transformation type.");
end
% 2-D transformation estimation
tformEst = imregcorr(moving_plane,Rmoving2D,fixed_plane,Rfixed2D,tf2D);
% modify the tformEst for skipping local minimum near zero when do
% translation, which comes from M-F transformation error
tl_core = tformEst.T(1:2, 1:2);
if isequal(reg_param.tform,"translation") && ...
        all(abs(eye(2)-tl_core)<1e-9,"all")
    tformEst.T(1:2, 1:2) = eye(2);
end

% !! first, fix z motion then adjust x-y motion
z_shift = 0;
% increase dimension of affine2d to affine3d
T = [[[tformEst.T(1:2,1:2),[0;0]];[0,0,1]];[tformEst.T(3,1:2),z_shift]];
T = [T,[0;0;0;1]];
tformEst = affine3d(T);

if nargout == 2
    % 1 memory copy from <imwarp>
    movolEst = imwarp(movol_denoised,Reg3D,tformEst,reg_param.interpRigid,...
        "OutputView",Reg3D,'SmoothEdge',true,'FillValues',bdv);
end

end

function plane = zProjection(vol,alg)
% maximum intensity as default

if ndims(vol) ~= 3
    error("invalid volume format: dimension not match");
end
if ~exist("alg","var")
    alg = 'max';
end
switch alg
    case 'avg'
        plane = mean(vol,3,"native");
    case 'max'
        plane = max(vol,[],3);
    case 'min'
        plane = min(vol,[],3);
    case 'std'
        plane = std(vol,0,3);
    case 'med'
        plane = median(vol,3);
    otherwise
        error("unsupported z-projection method");
end
end

% when do coarse registration, we need more robust for better perfermance
% the peak value of memory alloc is 2 fold of vol
function [vol,denoise_dim] = denoise3D(vol,varargin)

if ndims(vol)~=3
    error("invalid volume dimention: not 3-D");
end
[img_height,~,~] = size(vol);
% This function do pre-process before alignment
filterSize = [1,1,1];   % odd number: 3 elements vector
boundary = [0.002,0.998];
denoise_dim = [];
var_explained = 25;
% percentage, which need to cover the principle joint intensity
% and position distribution information

m = 1;
while m < numel(varargin)
    switch varargin{m}
        case 'filterSize'
            filterSize = varargin{m+1};
        case 'bound'
            boundary = varargin{m+1};
        case 'dd'
            denoise_dim = varargin{m+1};
        case 'CRC'
            var_explained = varargin{m+1};
    end
    m = m + 2;
end

if ~all(mod(filterSize,2)==1)
    error("Invalid filter size, must be odd.");
end
if numel(boundary)~=2 || ~isnumeric(boundary) || boundary(1)<0 ...
        || boundary(2)>1 || boundary(1)>=boundary(2)
    error("Invalid boundary, bad format.");
end

% 1. use quantile remove impulse noise
T_low = quantile(vol,boundary(1),"all");
T_high = quantile(vol,boundary(2),"all");
vol(vol<T_low)=T_low;
vol(vol>T_high)=T_high;

% 2. use filter denoise
% use median filter and box-filter
vol = medfilt3(vol,filterSize,'zeros');
vol = imboxfilt3(vol,filterSize,'padding',0);

% 3. use PCA extract basic structure
warning('off', 'stats:pca:ColRankDefX');    %give pca some drugs
% 'straighten' and 'folden' are inverse operator to each other
% straighten operator on 3-D tensor -> 2-D tensor:
optr_st = @(x)reshape(permute(x,[1,3,2]),[],size(x,2));
% folden operator on 2-D tensor -> 3-D tensor:
optr_fd = @(x,y)permute(reshape(x,y,size(x,1)/y,size(x,2)),[1,3,2]);

% change vol type as single for pca
% 1 fold memory alloc increase
vol = single(vol);
% permute and reconstruct array stand by...
vol = optr_st(vol);
if isempty(denoise_dim)
    % 'svd' for accuracy decompose
    [pc,score,~,~,explained,mu] = pca(vol,"Algorithm","svd");
    % select a variance explained threshold for better intensity estimate
    % and smaller memory loading
    denoise_dim = find(cumsum(explained)>=var_explained,1);
else
    % 'eig' for fast decompose
    [pc,score,~,~,~,mu] = pca(vol,"Algorithm","eig","NumComponents",denoise_dim);
end
vol = score(:,1:denoise_dim)*pc(:,1:denoise_dim)'+mu;

% 1 fold memory alloc decrease (recover)
vol = uint16(vol);

clearvars pc score mu;

botval = min(vol,[],"all");
% pixel value translation reconstruction
if botval<0
    vol = vol+abs(botval);
end
vol = optr_fd(vol,img_height);
end
