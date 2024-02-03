function [signal, tform] = reg_oc_auto(varargin)
%REGISTER3D_OC_AUTO This function help to align 3D volume series
%   This function can automatically align one channel volumetric images
%   NOTE: will be removed
%
%   reg = register3D_oc_auto()
%   [signal, tform] = register3D_oc_auto()
%   [signal, tform] = register3D_oc_auto(filename)
%   [signal, tform] = register3D_oc_auto(movinfo,filename)
%   [signal, tform] = register3D_oc_auto(___,Name,Value)
%
%   input:
%   - (optional) movinfo:  the structure of movie data(4/5-D) and
%                          information structure with field {'mov','opts'}
%   - (optional) filename: the name of file need to be aligned, only *.ims,
%                           *.nd2 and *.tiff is valid, image dimension: x*y*(c)*z*t
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
%       - RL: the region label, which could be
%           ["ORN","PN","LN","KC","MBON"], "ORN" as default
%       - GA: The pre-process transform parameter - gamma value, 2.0 as default
%       - InsTh: the intensity threshold for foreground detecting, 2.0 as default
%       - ScaleTh: the object size threshold (um) for foreground detecting, 3.0 as default
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
%       - GridUnit: the grid unit size when local registration, 
%                   ["auto", "1 1 1","2 2 1","3 3 1","4 4 1"], "auto" as default
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
%
%   see also: imregdemons, imregconfig, imregtform, imref3d, imregcorr,
%             pca, medfilt3, imboxfilt3, imhistmatchn, affine3d

% REGISTER3D_OC_AUTO:
% Version: 1.0.0
%   *** Basic registration functions
% Version: 1.1.0
%   *** append task manager for multi-user registration

% Copyright (c) 2022-2023, Weihan Li

% Validate input arguments
VALID_INTERP_METHOD = ["nearest", "linear", "cubic", "bilinear", "bicubic"];
VALID_REGMODE_METHOD = ["global-only", "local-only","mix"];
VALID_MODAL = ["multimodal", "monomodal"];
VALID_TFORM_MODE = ["translation", "rigid", "affine", "similarity"];
VALID_REGION_LABEL = ["ORN","PN","LN","KC","MBON"];
VALID_GRID_UNIT = ["auto","1 1 1","2 2 1","3 3 1","4 4 1"];

p = inputParser;
valid_file = @(x) isempty(x)||((isstring(x) || ischar(x)) && exist(x,"file"));
valid_movinfo = @(x) isstruct(x) && all(ismember(["mptr","opts"],string(fieldnames(x))));
valid_refvol = @(x)isstruct(x) && all(ismember(string(fieldnames(x)),["G","L"]));
valid_regframes = @(x) validateattributes(x,{'numeric'},{'row','positive','integer','increasing'});
valid_regmode = @(x)(ismember(x,VALID_REGMODE_METHOD));
valid_modal = @(x)(ismember(x,VALID_MODAL));
valid_tform = @(x)(ismember(x,VALID_TFORM_MODE));
valid_regionlabel = @(x)(ismember(x, VALID_REGION_LABEL));
valid_gridunit = @(x)(ismember(x, VALID_GRID_UNIT));
valid_initstep = @(x) strcmp(x,"auto") || (isscalar(x) && isreal(x) && x>0);
valid_maxitern = @(x)((isvector(x) && isPositiveIntegerValuedNumeric(x))...
    || strcmp(x,"auto"));
valid_minstep = @(x)(isscalar(x) && isreal(x) && x>0);
valid_itercoeff = @(x)validateattributes(x,{'numeric'},{'scalar','real','>',0});
valid_afs = @(x)validateattributes(x,{'numeric'},{'scalar','real','>=',0.5,'<=',3});
valid_gamma = @(x)validateattributes(x,{'numeric'},{'scalar','real','>=',0,'<=',4});
valid_insth = @(x)validateattributes(x,{'numeric'},{'scalar','real','>=',-5,'<=',5});
valid_scaleth = @(x)validateattributes(x,{'numeric'},{'scalar','real','>=',0});
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
default_initstep = "auto";           % the initial step when optimize registration
default_minstep = 1e-5;              % the minimum step when optimize registration
default_itercoeff = 0.5;             % the iteration coefficient
default_afs = 1.0;                   % the afs
default_max_iter_num = "auto";       % the maximum iteration number when optimize registration
default_region_label = "ORN";        % the default region label is "ORN"
default_grid_unit = "auto";          % the default grid unit for local registration
default_gamma = 2.0;                 % the default gamma value for contrast enhance
default_instah = 2.0;                % the intensity threshold for split foreground
default_scaleth = 2.0;               % the object size (um) low boundary, z-score
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
addParameter(p,'InitStep',default_initstep,valid_initstep);
addParameter(p,'MinStep',default_minstep,valid_minstep);
addParameter(p,'IterCoeff',default_itercoeff,valid_itercoeff);
addParameter(p,"RL",default_region_label, valid_regionlabel);
addParameter(p, "GA",default_gamma, valid_gamma);
addParameter(p, "InsTh", default_instah, valid_insth);
addParameter(p, "ScaleTh", default_scaleth, valid_scaleth);
addParameter(p, "GridUnit", default_grid_unit, valid_gridunit);
addParameter(p,'AFS',default_afs,valid_afs);
addParameter(p,'MaxIterN_Rigid',default_max_iter_num,valid_maxitern);
addParameter(p,'MaxIterN_NonRigid',default_max_iter_num,valid_maxitern);
addParameter(p,'VPL_Rigid',default_vpl,valid_vpl);
addParameter(p,'VPL_NonRigid',default_vpl,valid_vpl);
addParameter(p,'Interp_Rigid',default_interp,valid_interp);
addParameter(p,'Interp_NonRigid',default_interp,valid_interp);
addParameter(p,'LRTDS',default_lrtds,valid_lrtds);
addParameter(p,'BigFile',default_bigfile,valid_bigfile);
parse(p,varargin{:});

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

% pre-processing data that was RegFrames mentioned
if isequal(parser_results.RegFrames,intmax("uint16"))
    reg_frames = 1:opts.frames;
else
    reg_frames = intersect(parser_results.RegFrames, 1:opts.frames);
end
mov(:,:,:,:,setdiff(1:opts.frames,reg_frames)) = []; % remove the kept data
opts.frames = numel(reg_frames); % adjust the frames
opts.images = opts.frames*opts.slices*opts.channels; % adjust the total images

reg_param = struct( 'regmode',      parser_results.RegMode,...
    'fixvol',       [],...
    'RefVol',       parser_results.RefVol,...
    'modal',        parser_results.Modal,...
    'tform',        parser_results.Tform,...
    'RL',           parser_results.RL,...
    'GA',           parser_results.GA,...
    'InsTh',        parser_results.InsTh,...
    'ScaleTh',      parser_results.ScaleTh,...
    'GridUnit',     parser_results.GridUnit,...
    'initStep',     initstep,...
    'minStep',      minstep,...
    'iterCoeff',    itercoeff,...
    'afs',          parser_results.AFS,...
    'maxIterNumRigid',   max_iter_num_rigid,...
    'maxIterNumNonRigid',max_iter_num_nonrigid,...
    'vPLRigid',     parser_results.VPL_Rigid,...
    'vPLNonRigid',  parser_results.VPL_NonRigid,...
    'interpRigid',  parser_results.Interp_Rigid,...
    'interpNonRigid',parser_results.Interp_NonRigid,...
    'LRTDS',        parser_results.LRTDS,...
    'bigfile',      parser_results.BigFile);

% clear vars for debug easy
clearvars -except opts mov reg_param;

if ~isempty(mov)
    [tform, signal] = reg(mov, opts, reg_param);
else
    signal = [];
    tform = [];
    return;
end

    function [tform, mov] = reg(mov, opts, reg_param)
        % input:
        %   - ms_ptr: the movie signal channel (or disk file pointer)
        %   - opts:   the options of movie file
        %   - reg_param: the registration paramaters
        % output:
        %   - tform:  the geometric transformation
        %   - ms_ptr: the aligned movie signal channel (or disk file pointer)

        tform = cell(opts.frames,2);

        rm = GetRunningMode();

        % turn off the imregcorr 'weakpeak' warning
        warning('off','images:imregcorr:weakPeakCorrelation');

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

                    MOV_SIGBKG = quantile(mov, 0.05, "all");

                    % images pre-processing to uint8
                    mov_u8 = mov_preprocess(mov, reg_param.GA);

                    % generate global value: fixed volume and data scale
                    fv = GenFixVolProfile(mov_u8, reg_param.RefVol, "1X1X1");

                    % generate the raw reference coordinate system
                    RMOV_3D = imref3d(size(fv.fixvol_global_s_ds));

                    % extract the image border value, border width is 3
                    border_width = 3;

                    switch reg_param.regmode
                        case 'global-only'
                            bkgs_ud = [fv.fixvol_global_s_ds(1:border_width, :, :);...
                                fv.fixvol_global_s_ds(end-border_width+1:end, :, :)];
                            bkgs_lr = [fv.fixvol_global_s_ds(:, 1:border_width, :),...
                                fv.fixvol_global_s_ds(:, end-border_width+1:end, :)];
                            bkgs = single([bkgs_ud(:); bkgs_lr(:)]);
                            bkg_mu = mean(bkgs);
                            bkg_std = std(bkgs);

                            % generate the fixed square block
                            bb_fixed = generate_one_block(fv.fixvol_global_s_ds, ...
                                round(reg_param.ScaleTh/sqrt(opts.xRes*opts.yRes)));
                            if isempty(bb_fixed)
                                error("register3D_oc_auto:BadTemplate","Please check your template volume.");
                            end

                            % generate the blocked mov and reference
                            % coordiante system
                            [movie_fixed_rs, MOV_REF_FIXED] = movcrop(...
                                fv.fixvol_global_s_ds, bb_fixed);

                            for k = 1:opts.frames
                                if getappdata(bar,'canceling')
                                    break;
                                end
                                
                                % enhance and smooth movie_sig_u8
                                mov_u8_tmpv = mov_u8(:,:,:,k);
                                mov_u8_tmpv = medfilt3(mov_u8_tmpv,[5,5,1]);
                                mov_u8_tmpv = imadjustn(mov_u8_tmpv);

                                % set the background to 0
                                mov_u8_tmpv(mov_u8_tmpv < bkg_mu+reg_param.InsTh*bkg_std) = 0;

                                [tform_rough, bb] = coreg3D(mov_u8_tmpv,...
                                    bb_fixed, reg_param, opts);

                                [movie_sig_raw_rs, MOV_REF] = movcrop(...
                                    squeeze(mov(:,:,1,:,k)), bb);

                                tform{k,1} = imregtform(movie_sig_raw_rs, MOV_REF,...
                                    movie_fixed_rs, MOV_REF_FIXED,reg_param.tform,optimizer,metric,...
                                    "InitialTransformation",tform_rough);

                                tempv = squeeze(mov(:,:,1,:,k));
                                mov(:,:,1,:,k) = imwarp(tempv, tform{k,1}, ...
                                    "OutputView",RMOV_3D,"FillValues",MOV_SIGBKG);

                                % transform from double to single and downsampling
                                % for decreasing memory allocation
                                tform{k,2} = DownsamplingDisplacement(tform{k,2}, reg_param.LRTDS);

                                waitbar(k/opts.frames,bar,"register processing "...
                                    +num2str(k/opts.frames*100,3)+"% ...");
                            end

                        case 'local-only'
                            bkgs_ud = [fv.fixvol_local_s_ds(1:border_width, :, :);...
                                fv.fixvol_local_s_ds(end-border_width+1:end, :, :)];
                            bkgs_lr = [fv.fixvol_local_s_ds(:, 1:border_width, :),...
                                fv.fixvol_local_s_ds(:, end-border_width+1:end, :)];
                            bkgs = single([bkgs_ud(:); bkgs_lr(:)]);
                            bkg_mu = mean(bkgs);
                            switch reg_param.GridUnit   % should be an odd number
                                case "auto"
                                    gf = round(fix(size(mov_u8, 1)/128)/2)*2+1;
                                otherwise
                                    gf = str2double(reg_param.GridUnit.extract(1));   
                            end

                            % padding the movie to integer fold of grid_factor
                            paddings = (fix(size(mov_u8,[1,2])/gf)+1)*gf ...
                                - size(mov_u8,[1,2]);
                            paddings = (paddings + 3*mod(paddings, 2))/2;

                            mov_u8_pad = padarray(mov_u8,[paddings,0,0],bkg_mu,"both");

                            mov_reg_pad_ds = zeros([size(mov_u8_pad,[1,2])/gf, size(mov_u8_pad,[3,4])],...
                                "uint8");

                            for k = 1:size(mov_reg_pad_ds, 4)
                                mov_reg_pad_ds(:,:,:,k) = ...
                                    imresize3(mov_u8_pad(:,:,:,k), 'Scale', [1/gf, 1/gf, 1]);
                            end

                            clear mov_u8_pad;

                            fixed_ds = GenFixVolProfile(mov_reg_pad_ds, reg_param.RefVol, "1X1X1");
                            fixed_ds = fixed_ds.fixvol_local_s_ds;

                            for k = 1:size(mov_reg_pad_ds, 4)
                                if getappdata(bar,'canceling')
                                    break;
                                end

                                tmp_vol = mov_reg_pad_ds(:,:,:,k);

                                % do imregdemons
                                [tform_demons, ~] = imregdemons(tmp_vol,...
                                    fixed_ds, "AccumulatedFieldSmoothing",1.5,...
                                    "PyramidLevels",3, "DisplayWaitbar",false);

                                % scale displacement field representation
                                tform_demons = recov_tform(tform_demons, gf);

                                % crop displacement field representation
                                tform_demons = tform_demons(paddings(1)+1:end-paddings(1), ...
                                    paddings(2)+1:end-paddings(2), :, :);

                                % apply displacement field
                                tempv = squeeze(mov(:,:,1,:,k));
                                mov(:,:,1,:,k) = imwarp(tempv, tform_demons, ...
                                    "cubic","FillValues",bkg_mu);

                                waitbar(k/opts.frames,bar,"register processing "...
                                    +num2str(k/opts.frames*100,3)+"% ...");
                            end
                        case 'mix'
                            % do nothing
                            fprintf("No mixing pipeline for single core, process is skipped.");
                        otherwise
                    end
                    
                    if k == opts.frames,msg = "registration succeed.";else, ...
                            msg = "registration canceled.";end
                    waitbar(k/opts.frames,bar,msg);
                    pause(1);
                    delete(bar);
                otherwise
            end
        else
            switch reg_param.regmode
                case 'global-only'
                    bar = parwaitbar(opts.frames,'Waitbar',true);
                    [tform, mov] = align_global_cpu(mov,tform);
                case 'local-only'
                    bar = parwaitbar(opts.frames,'Waitbar',true);
                    if isequal(rm,'gpu') || isequal(rm,'multi-gpu')
                        [tform, mov] = align_local_gpu(mov,tform);
                    else
                        [tform, mov] = align_local_cpu(mov,tform);
                    end
                case 'mix'
                    bar = parwaitbar(2*opts.frames,'Waitbar',true);

                    [tform,mov] = align_global_cpu(mov,tform);

                    if isequal(rm,'gpu') || isequal(rm,'multi-gpu')
                        [tform,mov] = align_local_gpu(mov,tform);
                    else
                        [tform,mov] = align_local_cpu(mov,tform);
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

        function [t,ms_new] = align_global_cpu(ms,t)
            MOV_SIGBKG_PAR = quantile(ms, 0.05, "all");

            % images pre-processing to uint8
            mov_u8_par = mov_preprocess(ms, reg_param.GA);

            % generate global value: fixed volume and data scale
            fv_par = GenFixVolProfile(mov_u8_par, reg_param.RefVol, "1X1X1");

            % generate the raw reference coordinate system
            RMOV_3D_PAR = imref3d(size(fv_par.fixvol_global_s_ds));

            % extract the image border value, border width is 3
            border_width_par = 3;
            bkgs_ud_par = [fv_par.fixvol_global_s_ds(1:border_width_par, :, :);...
                fv_par.fixvol_global_s_ds(end-border_width_par+1:end, :, :)];
            bkgs_lr_par = [fv_par.fixvol_global_s_ds(:, 1:border_width_par, :),...
                fv_par.fixvol_global_s_ds(:, end-border_width_par+1:end, :)];
            bkgs_par = single([bkgs_ud_par(:); bkgs_lr_par(:)]);
            bkg_mu_par = mean(bkgs_par);
            bkg_std_par = std(bkgs_par);

            % generate the fixed square block
            bb_fixed_par = generate_one_block(fv_par.fixvol_global_s_ds, ...
                round(reg_param.ScaleTh/sqrt(opts.xRes*opts.yRes)));
            if isempty(bb_fixed_par)
                error("register3D_oc_auto:BadTemplate","Please check your template volume.");
            end

            % generate the blocked mov and reference
            % coordiante system
            [mov_fixed_rs_par, MOV_REF_FIXED_PAR] = movcrop(...
                fv_par.fixvol_global_s_ds, bb_fixed_par);

            task_mgr = TaskManager(opts, reg_param, opts.frames);
            [~, nw] = task_mgr.Allocate("cpu");

            parobj = OpenParpool(nw);
            
            parfor m = 1:opts.frames
                ms_m = squeeze(ms(:,:,1,:,m));
                tmpv = mov_u8_par(:,:,:,m);
                tmpv = medfilt3(tmpv,[5,5,1]);
                tmpv = imadjustn(tmpv);
                tmpv(tmpv < bkg_mu_par+reg_param.InsTh*bkg_std_par) = 0;

                % Set up the Initial Registration
                [optimizer_par,metric_par] = imregconfig(reg_param.modal);
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

                [tform_rough_par, bb_par] = coreg3D(tmpv,...
                    bb_fixed_par, reg_param, opts);

                [mov_rs_par, MOV_REF_PAR] = movcrop(ms_m, bb_par);

                tformm1 = imregtform(mov_rs_par, MOV_REF_PAR,...
                    mov_fixed_rs_par, MOV_REF_FIXED_PAR,reg_param.tform,...
                    optimizer_par,metric_par,"InitialTransformation",tform_rough_par);

                tempv_par = squeeze(ms(:,:,1,:,m));
                ms(:,:,1,:,m) = imwarp(tempv_par, tformm1, ...
                    "OutputView",RMOV_3D_PAR,"FillValues",MOV_SIGBKG_PAR);

                t(m,1) = {tformm1};
                bar.Send;
            end

            ms_new = ms;

            CloseParpool(parobj);

            task_mgr.Free();
        end

        function [t,ms_new] = align_local_gpu(ms,t)
            task_mgr = TaskManager(opts, reg_param, opts.frames);
            [~, nw] = task_mgr.Allocate("gpu");

            parobj = OpenParpool(nw);

            % generate the support variables
            [MOV_SIGBKG_PAR, mov_reg_pad_ds_par, fixed_par, gf_par, paddings_par]...
                = GenSupportVars(ms);

            paddings_ga = gpuArray(paddings_par);
            parfor m = 1:size(mov_reg_pad_ds_par, 4)
                fixed_ga = gpuArray(fixed_par);

                movie_ga = gpuArray(mov_reg_pad_ds_par(:,:,:,m));

                % do imregdemons
                [tform_demons_ga, ~] = imregdemons(movie_ga,...
                    fixed_ga, "AccumulatedFieldSmoothing",reg_param.afs,...
                    "PyramidLevels",reg_param.vPLNonRigid, "DisplayWaitbar",false); %#ok<*PFBNS>

                % scale displacement field representation
                tform_demons_ga = recov_tform(tform_demons_ga, gf_par);

                % crop displacement field representation
                tform_demons_ga = tform_demons_ga(paddings_ga(1)+1:end-paddings_ga(1), ...
                    paddings_ga(2)+1:end-paddings_ga(2), :, :);

                % apply displacement field
                tform_demons_par = gather(tform_demons_ga); % gather to memory
                tempv_par = squeeze(ms(:,:,1,:,m));
                ms(:,:,1,:,m) = imwarp(tempv_par, tform_demons_par, ...
                    "cubic","FillValues",MOV_SIGBKG_PAR);

                t(m,2) = {tform_demons_par};
                bar.Send;
            end

            for m = 1:opts.frames
                t{m, 2} = DownsamplingDisplacement(t{m, 2}, reg_param.LRTDS);
            end

            ms_new = ms;

            CloseParpool(parobj);

            task_mgr.Free();
        end

        function [t,ms_new] = align_local_cpu(ms,t)
            task_mgr = TaskManager(opts, reg_param, opts.frames);
            [~, nw] = task_mgr.Allocate("cpu");

            parobj = OpenParpool(nw);

            % generate the support variables
            [MOV_SIGBKG_PAR, mov_reg_pad_ds_par, fixed_ds_par, gf_par, paddings_par]...
                = GenSupportVars(ms);

            parfor m = 1:size(mov_reg_pad_ds_par, 4)
                tmp_vol_par = mov_reg_pad_ds_par(:,:,:,m);

                % do imregdemons
                [tform_demons_par, ~] = imregdemons(tmp_vol_par,...
                    fixed_ds_par, "AccumulatedFieldSmoothing",reg_param.afs,...
                    "PyramidLevels",reg_param.vPLNonRigid, "DisplayWaitbar",false); %#ok<*PFBNS>

                % scale displacement field representation
                tform_demons_par = recov_tform(tform_demons_par, gf_par);

                % crop displacement field representation
                tform_demons_par = tform_demons_par(paddings_par(1)+1:end-paddings_par(1), ...
                    paddings_par(2)+1:end-paddings_par(2), :, :);

                % apply displacement field
                tempv_par = squeeze(ms(:,:,1,:,m));
                ms(:,:,1,:,m) = imwarp(tempv_par, tform_demons_par, ...
                    "cubic","FillValues",MOV_SIGBKG_PAR);

                t(m,2) = {tform_demons_par};
                bar.Send;
            end

            for m = 1:opts.frames
                t{m, 2} = DownsamplingDisplacement(t{m, 2}, reg_param.LRTDS);
            end

            ms_new = ms;

            CloseParpool(parobj);

            task_mgr.Free();
        end

    end

    function [msp, mrpdp, fdp, gp, pp] = GenSupportVars(ms)
        msp = quantile(ms, 0.05, "all");

        % images pre-processing to uint8
        mov_u8_par = mov_preprocess(ms, reg_param.GA);

        % generate global value: fixed volume and data scale
        fv_par = GenFixVolProfile(mov_u8_par, reg_param.RefVol, "1X1X1");

        border_width_par = 3;
        bkgs_ud_par = [fv_par.fixvol_local_s_ds(1:border_width_par, :, :);...
            fv_par.fixvol_local_s_ds(end-border_width_par+1:end, :, :)];
        bkgs_lr_par = [fv_par.fixvol_local_s_ds(:, 1:border_width_par, :),...
            fv_par.fixvol_local_s_ds(:, end-border_width_par+1:end, :)];
        bkgs_par = single([bkgs_ud_par(:); bkgs_lr_par(:)]);
        bkg_mu_par = mean(bkgs_par);
        switch reg_param.GridUnit   % should be an odd number
            case "auto"
                gp = round(fix(size(mov_u8_par, 1)/128)/2)*2+1;
            otherwise
                gp = str2double(reg_param.GridUnit.extract(1));
        end

        % padding the movie to integer fold of grid_factor
        pp = (fix(size(mov_u8_par,[1,2])/gp)+1)*gp ...
            - size(mov_u8_par,[1,2]);
        pp = (pp + 3*mod(pp, 2))/2;

        mov_u8_pad_par = padarray(mov_u8_par,[pp,0,0],bkg_mu_par,"both");

        mrpdp = zeros([size(mov_u8_pad_par,[1,2])/gp, size(mov_u8_pad_par,[3,4])],...
            "uint8");

        for m = 1:size(mrpdp, 4)
            mrpdp(:,:,:,m) = ...
                imresize3(mov_u8_pad_par(:,:,:,m), 'Scale', [1/gp, 1/gp, 1]);
        end

        clear mov_u8_pad_par;

        fdp = GenFixVolProfile(mrpdp, reg_param.RefVol, "1X1X1");
        fdp = fdp.fixvol_local_s_ds;
    end

    function fv = GenFixVolProfile(mov, refvol, ds)
        mov = squeeze(mov);
        fixvol_global_s = GetFixedVol(mov, refvol.G);
        fixvol_local_s = GetFixedVol(mov, refvol.L);

        % get the image border fill value =========================
        borderval_global_s = imborderval(fixvol_global_s, [5,2,5,2,5,5]);
        borderval_global_s = mean(borderval_global_s([1,3,5,6]));
        borderval_local_s = imborderval(fixvol_local_s, [5,2,5,2,5,5]);
        borderval_local_s = mean(borderval_local_s([1,3,5,6]));

        fv.borderval_global_s = borderval_global_s;
        fv.borderval_local_s = borderval_local_s;

        switch ds
            case "auto"
                [fv.fixvol_global_s_ds, fv.ds_scale] = ReSample(fixvol_global_s);
                [fv.fixvol_local_s_ds, ~] = ReSample(fixvol_local_s, fv.ds_scale);
            otherwise
                [fv.fixvol_global_s_ds, fv.ds_scale] = ReSample(...
                    fixvol_global_s, 1/str2double(ds.extract(1)));
                [fv.fixvol_local_s_ds, ~] = ReSample(...
                    fixvol_local_s, fv.ds_scale);
        end

        % balance and enhance
        fv.fixvol_global_s_ds = medfilt3(fv.fixvol_global_s_ds, [5,5,1]);
        fv.fixvol_global_s_ds = imadjustn(fv.fixvol_global_s_ds);
        fv.fixvol_local_s_ds = medfilt3(fv.fixvol_local_s_ds, [5,5,1]);
        fv.fixvol_local_s_ds = imadjustn(fv.fixvol_local_s_ds);

    end
end

function tform = recov_tform(tform_ds, sf, alg)
arguments
    tform_ds (:,:,:,:) double;
    sf (1,1) {mustBeReal, mustBePositive, mustBeFinite} = 1;
    alg (1,1) string {mustBeMember(alg, ["linear","nearest","cubic", ...
        "makima","spline"])} = "linear";
end

tform = zeros([size(tform_ds, [1,2])*sf, size(tform_ds, [3,4])],"like",tform_ds);

dispF_X = tform_ds(:,:,:,1);
dispF_Y = tform_ds(:,:,:,2);
dispF_Z = tform_ds(:,:,:,3);

[X, Y, Z] = meshgrid(sf-1:sf:size(tform,2)-1, ...
                     sf-1:sf:size(tform,1)-1, ...
                     1:size(tform,3));

[Xq, Yq, Zq] = meshgrid(1:size(tform,2), 1:size(tform,1), 1:size(tform, 3));

tform(:,:,:,1) = interp3(X, Y, Z, dispF_X, Xq, Yq, Zq, alg);
tform(:,:,:,2) = interp3(X, Y, Z, dispF_Y, Xq, Yq, Zq, alg);
tform(:,:,:,3) = interp3(X, Y, Z, dispF_Z, Xq, Yq, Zq, alg);

% refine the border
tform(isnan(tform)) = 0;
end


function rect_info = generate_one_block(V, r, w)
J = imerode(V,strel("cube",3));

% dilate enhance connected inner and clear border
J = imdilate(medfilt3(imclearborder(J), [5,5,5]), strel("square",5));

% generate BW volume
BW = bwareaopen(imbinarize(J,"global"),round(pi*4/3*r^3),26);

% append border on BW
BW = padarray(BW,[10,10,0],0,"both");

% extract the structure with optimized value
stats = regionprops3(BW,"Volume","SurfaceArea","Solidity","BoundingBox");

if isempty(stats)
    warning("register3D_oc_auto:ObjectLost","Please check the pre-processing filter.");
    rect_info = [];
    return;
end

% select the maximum volume and Spherity
% or we trained a network for neural regions detection
spherity = ((4*pi)^(1/3)*(3*stats.Volume).^(2/3))./stats.SurfaceArea;
[~, rp] = max(stats.Volume.*spherity.*stats.Solidity);
stats = stats(rp, :);

% generate masking squared
if ~exist("w","var")
    [~, pos] = min(stats.BoundingBox(4:5));
    stats.BoundingBox(pos) = stats.BoundingBox(pos) ...
        - abs(diff(stats.BoundingBox(4:5)))/2;
    stats.BoundingBox(pos+3) = stats.BoundingBox(pos+3) ...
        + abs(diff(stats.BoundingBox(4:5)))/2;
else
    % constrain the width and height with w
    % expand or shrink the boundary
    stats.BoundingBox(1) = stats.BoundingBox(1) ...
        + (w - stats.BoundingBox(4))/2;
    stats.BoundingBox(2) = stats.BoundingBox(2) ...
        + (w - stats.BoundingBox(5))/2;
    stats.BoundingBox(4) = w;
    stats.BoundingBox(5) = w;
end

rect_info = round([stats.BoundingBox([1,2]), ...
                       stats.BoundingBox(3), ...
                       stats.BoundingBox(4), ...
                       stats.BoundingBox(6)]);
end

% About maximum 3 fold memory of movol used <coreg3D>
function [tformEst,bb] = coreg3D(movol,bb_fixed,reg_param,opts)
% this function register movol to fixvol by using morphology and 
% object-matching based on ANN
% Input:
%   - movol: 3D marix
%   - bb_fixed: 1 X 5 fixed volume square border
%   - ...

bb = generate_one_block(movol,...
    round(reg_param.ScaleTh/sqrt(opts.xRes*opts.yRes)),bb_fixed(4));
if isempty(bb)
    warning("register3D_oc_auto:BadMotion", "Registration skipped.");
end

% construct translation matrix and appling 3D shift
x_shift = bb_fixed(1) - bb(1);
y_shift = bb_fixed(2) - bb(2);
z_shift = (bb_fixed(3)+bb_fixed(5)/2)- (bb(3)+bb(5)/2);
A_3D = [[eye(3),[0;0;0]];[x_shift,y_shift,z_shift,1]];
tformEst = affine3d(A_3D);
end

function [MOVRS, MOVREF] = movcrop(mov, bb)
movie_fixed_ROI = mov(bb(2):bb(2)+bb(4)-1, bb(1):bb(1)+bb(4)-1,:);
MOVRS = imresize3(movie_fixed_ROI,[128,128,size(movie_fixed_ROI,3)]);
MOVREF = imref3d(size(MOVRS), ...
    [bb(1)-0.5,bb(1)+bb(4)-1.5], ...
    [bb(2)-0.5,bb(2)+bb(4)-1.5] ,...
    [0.5, size(MOVRS,3)+0.5]);
end

function mov_u8 = mov_preprocess(mov_u16, ga)
% images pre-processing
movie_sig = squeeze(mov_u16);

% shrink LUT range for robustness imaging
intensity_upb = quantile(movie_sig,0.999,"all");
movie_sig(movie_sig>intensity_upb) = intensity_upb;
intensity_lob = min(movie_sig,[],"all");

mov_u8 = uint8(single(movie_sig - intensity_lob)./single(intensity_upb-intensity_lob)*255);

clearvars movie_sig intensity_lob intensity_upb;

% using gamma transformation for contrast enhance
mov_u8 = uint8((double(mov_u8)/255).^ga*255);
end