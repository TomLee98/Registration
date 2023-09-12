function [signal, tform, marker] = register3D_lctc_auto(varargin)
%REGISTER3D_LCTC This function (long-term chain with two channels) uses point 
% clouds and registration chain do long term 3d registration, note that 
% the registration mode uses rigid、affine and non-rigid algorithms on the whole image
%   Note that: consider the long-term experiment, we always need RTIO
%   (Real time input output) method but waste some time
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
%       - KeyFrames: the key frames series, string with points or "auto"
%       - IntensityThreshold: the lower intensity threshold for
%                           pre-foreground extraction, 1~100, 97 as default
%       - ScaleThreshold: the scale threshold for pre-foreground
%                        extraction, omit little objects for robust registration
%       - DS_PointCloud: the downsampling options for point cloud
%                      registration, could be "GA","RAND","NU", "GA" as default
%       - DS_PointCloud_Param: the downsampling parameter value
%       - OutlierRatio: the outlier ratio for more anti-noise in points 
%                       cloud registration, 0~1, 0.1 as default
%       - MaxIterN_PointCloud: the maximum iteration number when points 
%                              cloud registration is running, which is one 
%                              of the stop condition, 50 as default
%       - ErrorLimit: the error limitation when points cloud is running,
%                     which is one of the stop condition
%       - InitStep: the initial step when optimizal two volumes loss, it is
%                   initial radius(default: 6.25e-3) when using "multimodal", or
%                   maximum step length(default: 6.25e-2) when using "monomodal"
%       - MinStep: the minimum iteration step when optimizal two volumes
%                  loss, which meature the convergance
%       - DS_Voxel: the downsampling options for speed up and more robust, could
%                   be "auto", "2X2", "3X3", "4X4", "auto" as default
%       - MaxIterN_Voxel: the maximum iteration number when alignment, 100 (default)
%       - IterCoeff: the iteration coefficient for loop, if Modal is monomodal,
%                    then the range of IterCoeff is (0,1), 0.5 as default;
%                    and if Modal is multimodal, then the range of
%                    IterCoeff is (0,+inf), 1.05 as default
%       - AutoContrast: the flag for auto contrast on marker channel, for
%                       better local registration, true(default)/false
%       - CompAcc: the compensate accuracy, the imhistmatchn binning counts
%       - AFS: the accumulative field smoothing, the AFS bigger and the
%              displacement is more smooth, usual 0.5-3, 1 as default
%       - Chla: the channel you want to align, this section is enabled
%                only if multi-channel ["r"(default),"g","b"]
%       - Chls: the channel you want to extract signal ["r","g"(default),"b"]
%       - LRTDS: the local registration transformation downsampling level,
%                integer number for how many pixel cause a sampling point,
%                1(no downsample) as default
%   output:
%   - tform: 1st-col: the global transform object, affine2d or affine3d
%            2nd-col: the local displacement field, double
%   - signal: the signal channel has been aligned, 4-D uint16
%   - marker: the marker channel has been aligned, 4-D uint16
%
%   see also: imregdemons, imregconfig, imregtform, imref3d, pcdownsample,
%             pcregistercpd, imwarp

% REGISTER3D_TC_AUTO:
% Version: 1.0.0
%   *** Basic registration functions

% Copyright (c) 2022-2023, Weihan Li

VALID_PC_DS_OPTIONS = ["GA","RAND","NU"];
VALID_VX_DS_OPTIONS = ["auto","1X1","2X2","3X3"];
VALID_CHANNEL_LABEL = ["r", "g", "b"];
DEFAULT_CHAIN_LENGTH = 20;

p = inputParser;
valid_file =        @(x)assert(isempty(x)||((isstring(x) || ischar(x)) && exist(x,"file")));
valid_movinfo =     @(x)assert(isstruct(x) && all(ismember(["mov","opts"],string(fieldnames(x)))));
valid_refvol =      @(x)assert(isstruct(x) && all(ismember(string(fieldnames(x)),["G","L"])));
valid_regframes =   @(x)validateattributes(x,{'numeric'},{'row','positive','integer','increasing'});
valid_keyframes =   @(x)validateattributes(x, {'string'},{'scalar'});
valid_ityth =       @(x)validateattributes(x, {'numeric'},{'scalar','>=',0,'<=',100});
valid_scath =       @(x)validateattributes(x, {'numeric'},{'scalar','>=',0.1,'<=',1000});
valid_pc_ds =       @(x)assert(ismember(x, VALID_PC_DS_OPTIONS));
valid_vx_ds =       @(x)assert(ismember(x, VALID_VX_DS_OPTIONS));
valid_ds_param =    @(x)validateattributes(x, {'numeric'}, {'scalar','positive'});
valid_olr_ratio =   @(x)validateattributes(x,{'numeric'},{'scalar','real','>=',0,'<',1});
valid_pc_maxitern = @(x)validateattributes(x,{'numeric'},{'scalar','integer','>=',1});
valid_err_limit =   @(x)validateattributes(x,{'numeric'},{'scalar','real','>',0});
valid_initstep =    @(x)assert(isscalar(x) && isreal(x) && x>0);
valid_minstep =     @(x)assert(isscalar(x) && isreal(x) && x>0);
valid_vx_maxitern = @(x)assert(isvector(x) && isPositiveIntegerValuedNumeric(x));
valid_itercoeff =   @(x)validateattributes(x,{'numeric'},{'scalar','real','>',0,'<',1});
valid_autocontrast= @(x)validateattributes(x,{'logical'},{'scalar'});
valid_compacc =     @(x)assert(ismember(x,VALID_COMP_ACC));
valid_afs =         @(x)validateattributes(x,{'numeric'},{'scalar','real','>=',0.5,'<=',3});
valid_chla =        @(x)assert(isscalar(x) && (isstring(x)||ischar(x)) ...
                                && ismember(lower(x),VALID_CHANNEL_LABEL));
valid_chls =        @(x)assert(isscalar(x) && (isstring(x)||ischar(x)) ...
                                && ismember(lower(x),VALID_CHANNEL_LABEL));

%======================== DEFAULT PARAMETER SETTING =======================
default_filename = [];
default_movinfo = struct('mov',[],'opts',[]);
default_refvol = struct("G",   ["mean", "1"], ...
                        "L",   ["mean", "1"]);
default_regframes = intmax("uint16");
default_keyframes = "auto";
default_ityth = 97;
default_scath = 3;
default_pc_ds = "GA";
default_vx_ds = "AUTO";
default_ds_param = 2;
default_olr_ratio = 0.1;
default_pc_maxitern = 50;
default_err_limit = 1e-5;
default_initstep = 6.25e-2;
default_minstep = 1e-5;
default_vx_maxitern = [100,50,25];
default_itercoeff = 0.5;
default_autocontrast = true;
default_compacc = 1024;
default_afs = 1.0;
default_channel_align = "r";
default_channel_signal = "g";

%==========================================================================

addOptional(p,'movinfo',default_movinfo,valid_movinfo);
addOptional(p,'filename',default_filename,valid_file);
addParameter(p,'RefVol',default_refvol,valid_refvol);
addParameter(p,'RegFrames',default_regframes,valid_regframes);
addParameter(p,'KeyFrames',default_keyframes,valid_keyframes);
addParameter(p,'IntensityThreshold',default_ityth,valid_ityth);
addParameter(p,'ScaleThreshold',default_scath, valid_scath);
addParameter(p,'DS_PointCloud',default_pc_ds,valid_pc_ds);
addParameter(p,'DS_Voxel',default_vx_ds,valid_vx_ds);
addParameter(p,'DS_PointCloud_Param',default_ds_param,valid_ds_param);
addParameter(p,'Outlier_Ratio',default_olr_ratio,valid_olr_ratio);
addParameter(p,'MaxIterN_PointCloud',default_pc_maxitern,valid_pc_maxitern);
addParameter(p,'ErrorLimit',default_err_limit,valid_err_limit);
addParameter(p,'InitStep',default_initstep,valid_initstep);
addParameter(p,'MinStep',default_minstep,valid_minstep);
addParameter(p,'MaxIterN_Voxel',default_vx_maxitern,valid_vx_maxitern);
addParameter(p,'IterCoeff',default_itercoeff,valid_itercoeff);
addParameter(p,'AutoContrast',default_autocontrast,valid_autocontrast);
addParameter(p,'CompAcc',default_compacc,valid_compacc);
addParameter(p,'AFS',default_afs,valid_afs);
addParameter(p,'Chla',default_channel_align,valid_chla);
addParameter(p,'Chls',default_channel_signal,valid_chls);
parse(p, varargin{:});

% because the different memory management between
% unix(struct with deep copy) and windows(struct with shallow copy)
pr = p.Results;

% load the file
if isempty(pr.movinfo.mov)
    [opts, mov,~,~] = loadfile(pr.filename);
else
    mov = pr.movinfo.mov;
    opts = pr.movinfo.opts;
end

% resetup some parameters
switch pr.DS_Voxel
    case "AUTO"
        ds_vx = 64./max(opts.width, opts.height);
    case "2X2"
        ds_vx = 0.5;
    case "3X3"
        ds_vx = 1/3;
    case "OFF"
        ds_vx = 1;
end

if pr.KeyFrames == "auto"
    % select the frames with same blank but shift half
    key_frs_idx = round(linspace(1, opts.frames, DEFAULT_CHAIN_LENGTH+1));
    key_frs_idx = key_frs_idx(1:end-1) + fix(diff(key_frs_idx(1:2)));
else
    key_frs = pr.KeyFrames.replace("end",num2str(opts.frames));
    key_frs_idx = str2num(key_frs); %#ok<ST2NM>
end

c_order = opts.cOrder;
% generate global value: fixed volume and data scale
[MA_MIN, MA_MAX] = getMinMaxIn(mov, pr.Chla, c_order);
fv = genFixVolProfile(mov, pr, c_order, ds_vx);

% extract the key frames
key_mov = squeeze(mov(:,:,pr.Chla==c_order,:,key_frs_idx));

% pre-processing data that was RegFrames mentioned
if isequal(pr.RegFrames,intmax("uint16"))
    reg_frames = 1:opts.frames;
else
    reg_frames = intersect(pr.RegFrames, 1:opts.frames);
end
mov(:,:,:,:,setdiff(1:opts.frames,reg_frames)) = []; % remove data
opts.frames = numel(reg_frames);                     % adjust the frames
opts.images = opts.frames*opts.slices*opts.channels; % adjust the total images

reg_param = struct('fixvol',       fv,...
                    'regFrames',    reg_frames,...
                    'keyFramesIdx', key_frs_idx,...
                    'intensityTh',  pr.IntensityThreshold,...
                    'scaleTh',      pr.ScaleThreshold,...
                    'dsPC',         pr.DS_PointCloud,...
                    'dsParam',      pr.DS_PointCloud_Param,...
                    'olrRatio',     pr.OutlierRatio,...
                    'maxIterPC',    pr.MaxItern_PointCloud,...
                    'errorLimit',   pr.ErrorLimit,...
                    'initStep',     pr.InitStep,...
                    'minStep',      pr.MinStep,...
                    'dsVX',         ds_vx,...
                    'maxIterVX',    pr.MaxIterN_Voxel,...
                    'iterCoeff',    pr.IterCoeff,...
                    'autoContrast', pr.AutoContrast,...
                    'compAcc',      pr.CompAcc,...
                    'afs',          pr.AFS,...
                    'chlA',         pr.Chla,...
                    'chlS',         pr.Chls,...
                    'chlMode',      c_order);

clearvars -except opts mov key_mov fv reg_param MA_MIN MA_MAX;

chl = [find(reg_param.chlA==reg_param.chlMode), ...
                   find(reg_param.chlS==reg_param.chlMode)];
[ma_ptr, ms_ptr] = GenFilePointer(mov, chl);
clearvars mov;

% generate the registration chain key point
[tf_affine, tf_nonrigid, ~] = gen_reg_kp(key_mov, fv.fixvol_global_a, opts);
clearvars fv;

if nargout == 1
    tform = reg(ma_ptr, ms_ptr, tf_affine, tf_nonrigid, opts, reg_param);
elseif nargout == 2
    [tform, signal] = reg(ma_ptr, ms_ptr, tf_affine, tf_nonrigid, opts, reg_param);
elseif nargout == 3
    [tform, signal, marker] = reg(ma_ptr, ms_ptr, tf_affine, tf_nonrigid, opts, reg_param);
end

    function [tform, ms_ptr, ma_ptr, good_flag] = reg(ma_ptr, ms_ptr, tf_affine, tf_nonrigid)
        %   (1) register to key point in the chain
        %   (2) apply the chain shift
        %   (3) evaluate the registration result, do twice registration on bad part
        n_out = nargout;

        [rm, wkn] = getRunningMode();

        if strcmp(rm, 'cpu')
            tform = cell.empty(0, 2);
            warning("register3D_lctc_auto:tooWeakComputationResource",...
                "Please use multicore CPU for long-term registration.");
            return;
        end

        tform = cell(opts.frames, 2);
        smse = nan(opts.frames, 1);
        good_flag = nan(opts.frames, 1);

        %  read data as more as possible for speed up
        block_n = getCpuBlockNumber(opts);
        load_loop_n = ceil(opts.frames/block_n);
        openParpool(min(block_n, wkn)); % rigid/affine with cpu first

        % extract the params for avoiding data broadcast
        reg_frame = reg_param.regFrames;
        key_frame_idx = reg_param.keyFramesIdx;
        pc_ds = reg_param.ds_PC;
        pc_ds_param = reg_param.dsParam;
        pc_olr = reg_param.olrRatio;
        pc_max_iter_n = reg_param.maxIterPC;
        vx_ds = reg_param.ds_VX;
        vx_rlx_factor = reg_param.iterCoeff;
        vx_max_iter_n = reg_param.maxIterVX;
        vx_max_step_len = reg_param.initStep;
        vx_min_step_len = reg_param.minStep;
        tol = reg_param.errorLimit;
        fine_tune_max_n = 3;
        rs = [opts.xRes, opts.yRes, opts.zRes];
        RA = imref3d([opts.height, opts.width, opts.slices], ...
            rs(1), rs(2), rs(3));
        RA_ds = imref3d([opts.height, opts.width, opts.slices], ...
            rs(1)/vx_ds, rs(2)/vx_ds, rs(3));

        fix_init = reg_param.fixvol.fixvol_local_a;
        fix_init_ds = reg_param.fixvol.fixvol_global_a_ds;
        std_const = prctile(fix_init, 99, "all");

        % cpu registration:
        for k = 1:load_loop_n
            block_frames = min(block_n, opts.frames-(k-1)*block_n);
            st_idx = (k-1)*block_n + 1;
            ed_idx = (k-1)*block_n + block_frames;

            % load the temp data to memory
            ma = ma_ptr(:,:,:,st_idx:ed_idx);
            ms = ms_ptr(:,:,:,st_idx:ed_idx);
            smse_block = nan(block_frames, 1);
            flag_block =  nan(block_frames, 1);
            tf_affine_block = cell(block_frames, 1);

            % using parallel computation for speed up
            parfor n = 1:block_frames
                [optimizer, metric] = imregconfig("monomodal");
                optimizer.MaximumIterations = vx_max_iter_n
                optimizer.MinimumStepLength = vx_min_step_len;
                optimizer.MaximumStepLength = vx_max_step_len;
                optimizer.RelaxationFactor = vx_rlx_factor;

                % select the key frame at first
                [~, pos] = min(abs(key_frame_idx - reg_frame(st_idx+n-1))); %#ok<*PFBNS>
                key_frame_affine = tf_affine{pos};

                % 1. rigid registration based on point cloud
                fixvol = key_mov(:,:,:,pos);
                movol_ma = ma(:,:,:,n);
                movol_ds = DownSampling(movol_ma, vx_ds);
                pts_fix = vol2pts(fixvol, rs);
                pts_mov = vol2pts(movol_ma, rs);

                % using cloud registration cpd -> affine registration
                % downsampling, rough registration, make sure the structured space
                switch pc_ds
                    case "GA"
                        pts_fix_ds = pcdownsample(pts_fix, "gridAverage", pc_ds_param);
                        pts_mov_ds = pcdownsample(pts_mov, "gridAverage", pc_ds_param);
                    case "RAND"
                        pts_fix_ds = pcdownsample(pts_fix, "random", pc_ds_param);
                        pts_mov_ds = pcdownsample(pts_mov, "random", pc_ds_param);
                    case "NU"
                        pts_fix_ds = pcdownsample(pts_fix, "nonuniformGridSample", pc_ds_param);
                        pts_mov_ds = pcdownsample(pts_mov, "nonuniformGridSample", pc_ds_param);
                    otherwise
                end

                % first rigid estimation
                tfs_affine = pcregistercpd(pts_mov_ds, pts_fix_ds, ...
                    "Transform","Rigid", "OutlierRatio",pc_olr,...
                    "MaxIterations",pc_max_iter_n, "Tolerance",tol);

                % second affine estimation
                tfs_affine = imregtform(movol_ds, RA_ds, fix_init_ds, RA_ds, "affine", ...
                    optimizer, metric, "InitialTransformation", tfs_affine);

                % using key frame affine as an estimation
                tfs_affine_fi = affinetform3d();
                tfs_affine_fi.A = tfs_affine.A*key_frame_affine.A;

                % apply the 'final' registration and evaluate the results
                MOV_BKG_MA = prctile(movol_ma, 10, "all");
                ma_aligned = imwarp(movol_ma, RA, tfs_affine_fi, "linear",...
                            "OutputView",RA,"FillValues",MOV_BKG_MA);
                smse_block(n) = immse(ma_aligned, fix_init)/std_const;

                % update the optimizer to avoid the possible bad convergence
                optimizer.MinimumStepLength = reg_param.minStep/10;
                optimizer.MaximumStepLength = reg_param.initStep/5;
                
                tune_times = 0;
                flag_block(n) = (smse_block(n) <= 0.5);
                % 0.5 as intensity iteration hard threshold
                while ~flag_block(n) && (tune_times < fine_tune_max_n)
                    % second update tforms for fine tuning
                    tfs_affine_fi = imregtform(movol_ma, RA, fix_init, RA, "affine", ...
                        optimizer, metric, "InitialTransformation", tfs_affine_fi);
                    ma_aligned = imwarp(movol_ma, RA, tfs_affine_fi, "linear",...
                        "OutputView",RA,"FillValues",MOV_BKG_MA);
                    smse_block(n) = immse(ma_aligned, fix_init)/std_const;
                    flag_block(n) = (smse_block(n) <= 0.5);
                    tune_times = tune_times + 1;
                end

                switch n_out
                    case 2
                        % apply affine on signal
                        movol_ms = ms(:,:,:,n);
                        MOV_BKG_MS = prctile(movol_ms, 10, "all");
                        ms(:,:,:,n) = imwarp(movol_ms, RA, tfs_affine_fi, "linear",...
                            "OutputView",RA,"FillValues",MOV_BKG_MS);
                    case 3
                        % apply affine on both signal and marker
                        movol_ms = ms(:,:,:,n);
                        MOV_BKG_MS = prctile(movol_ms, 10, "all");
                        ms(:,:,:,n) = imwarp(movol_ms, RA, tfs_affine_fi, "linear",...
                            "OutputView",RA,"FillValues",MOV_BKG_MS);
                        ma(:,:,:,n) = ma_aligned;
                    otherwise
                end

                tf_affine_block{n} = tfs_affine_fi;
            end

            tform(st_idx:ed_idx, 1) = tf_affine_block;
            smse(st_idx:ed_idx) = smse_block;
            good_flag(st_idx:ed_idx) = flag_block;

            % modify the hard drive data
            switch n_out
                case 2
                    ms_ptr(:,:,:,st_idx:ed_idx) = ms;
                case 3
                    ms_ptr(:,:,:,st_idx:ed_idx) = ms;
                    ma_ptr(:,:,:,st_idx:ed_idx) = ma;
                otherwise
            end
        end

        if strcmp(rm, 'multi-cpu')
            [tform_nrd, ms_ptr, ma_ptr] = imregdemons_fast_cpu(ma_ptr, ms_ptr);
        else
            [tform_nrd, ms_ptr, ma_ptr] = imregdemons_fast_gpu(ma_ptr, ms_ptr);
        end

        tform(:, 2) = tform_nrd;

    end

    function [tform, ms_ptr, ma_ptr] = imregdemons_fast_cpu(ma_ptr, ms_ptr)

    end

    function [tform, ms_ptr, ma_ptr] = imregdemons_fast_gpu(ma_ptr, ms_ptr)
        
    end


    % TODO: the tf_affine compatibility problem
    function [tf_affine, tf_nonrigid, rmses, exit_flag] = gen_reg_kp(mov, fix_init, opts)
        % This function generate the key point of registration chain
        % each point will be registered to the settle template
        % Note that: the error propagation will be significant when the
        % chain is too long
        % Input:
        %   - mov: 4-D uint16 array, the move volume series, [x,y,z,t]
        %   - fix_init: 3-D uint16 array, the fixed volume, [x,y,z]
        %   - opts: 1-by-8 table, channels and frames are wrong because crop
        % Output:
        %   - tf_affine: t-by-1 cell array, element is affinetform3d or
        %       affine3d object (MATLAB < 9.13, TODO)
        %   - tf_nonrigid: t-by-1 cell array, element is [x,y,z,3]
        %       displacement field array
        %   - rmses: t-by-1 array, point cloud registration rmse
        %   - exit_flag: logical, true for normal exit, false for breaking

        % expend the mov
        mov = cat(4, fix_init, mov);
        vol_size = size(mov, 1:3);
        frames = size(mov, 4);

        % setup the output vars
        tf_affine = cell(frames, 1);
        tf_affine{1} = affinetform3d(eye(4));
        tf_nonrigid = cell(frames, 1);
        tf_nonrigid{1} = zeros([vol_size, 3]);
        rmses = zeros(frames, 1);
        exit_flag = true;

        % setup temporary vars
        res_arr = [opts.xRes, opts.yRes, opts.zRes];
        RA = imref3d(vol_size, res_arr(1), res_arr(2), res_arr(3));
        fix_init_ds = DownSampling(fix_init, 0.25);
        RA_ds = imref3d(size(fix_init_ds), res_arr(1)/0.25, res_arr(2)/0.25, res_arr(3));
        [optimizer, metric] = imregconfig("monomodal");
        optimizer.MinimumStepLength = reg_param.minStep;
        optimizer.MaximumStepLength = reg_param.initStep;
        optimizer.RelaxationFactor = reg_param.iterCoeff;

        for k = 1:frames-1
            % generate points cloud
            fixvol = mov(:,:,:,k);
            movol = mov(:,:,:,k+1);

            pts_fix = vol2pts(fixvol, res_arr);
            pts_mov = vol2pts(movol, res_arr);
            if pts_fix.Count <= 5 || pts_mov.Count <= 5
                exit_flag = false;
                break;
            end

            movol_ds = DownSampling(movol, reg_param.ds_vx);
            MOV_BKG = prctile(movol, 10, "all");

            % using cloud registration cpd -> affine registration
            % downsampling, rough registration, make sure the structured space
            switch reg_params.dsPC
                case "GA"
                    pts_fix_ds = pcdownsample(pts_fix, "gridAverage", reg_param.dsParam);
                    pts_mov_ds = pcdownsample(pts_mov, "gridAverage", reg_param.dsParam);
                case "RAND"
                    pts_fix_ds = pcdownsample(pts_fix, "random", reg_param.dsParam);
                    pts_mov_ds = pcdownsample(pts_mov, "random", reg_param.dsParam);
                case "NU"
                    pts_fix_ds = pcdownsample(pts_fix, "nonuniformGridSample", reg_param.dsParam);
                    pts_mov_ds = pcdownsample(pts_mov, "nonuniformGridSample", reg_param.dsParam);
                otherwise
            end

            % first estimation for key k -> key k+1 points cloud registration
            [tform_affine, ~, rmses(k+1)] = pcregistercpd(pts_mov_ds, pts_fix_ds, ...
                "Transform","Rigid","OutlierRatio",0.1,"MaxIterations",100,"Tolerance",1e-7);
            tf_affine{k+1} = affinetform3d();

            % update the fixed -> k+1 transformation estimation
            tf_affine{k+1}.A = tform_affine.A*tf_affine{k}.A;

            % second update tforms again to avoid accumulation of errors
            tf_affine{k+1} = imregtform(movol_ds, RA_ds, fix_init_ds, RA_ds, "affine", ...
                optimizer, metric, "InitialTransformation", tf_affine{k+1});

            % update the fixed -> k+1 fine tuned nonrigid estimation
            reg_affined = imwarp(movol, RA, tf_affine{k+1},"linear",...
                "OutputView",RA,"FillValues",MOV_BKG);
            [tf_nonrigid{k+1}, ~] = imregdemons(reg_affined, fix_init, [100,50,25],...
                "AccumulatedFieldSmoothing",1,"DisplayWaitbar",false);
        end

        % remove the first(init) transformation
        tf_affine(1) = [];
        tf_nonrigid(1) = [];
    end



end




function fv = genFixVolProfile(mov, ps, cd, ds)
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

[fv.fixvol_global_a_ds, fv.ds_scale] = DownSampling(fixvol_global_a, ds);
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

function parobj = openParpool(n)
% generate parcluster
% modify
pcl = parcluster("Reg3D_Server");
% make sure your cpu supports 'hyper-threads' technology
pcl.NumThreads = 2;
% get present parpool
parobj = gcp("nocreate");

if isempty(parobj)
    parobj = parpool(pcl, [1,n], 'SpmdEnabled',false);
elseif parobj.NumWorkers ~= n
    % restart parpool
    delete(gcp);
    parobj = parpool(pcl, [1,n], 'SpmdEnabled',false);
end
end

function closeParpool(parobj)
% close parpool
delete(parobj);
end

function gn = getGpuBlockNumber(opts)
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

function cn = getCpuBlockNumber(opts)
% This function get the avaiable memory size and give a bench
% size of a processing job, which is the frame could be loaded
UINT16_BYTES = 2;
FOLD_RATIO = 5;
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
end
if cn < 1
    error('No available enough memory.');
end
end