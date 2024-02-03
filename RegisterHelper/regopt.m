classdef regopt
    %REGOPT This class is registration option defination, which is packaged
    %value class
    
    properties(Access=private, Hidden)
        % ============= one/two channel(s) common properties ==============
        reg_modal       (1,1) string {mustBeMember(reg_modal, ["multimodal", "monomodal"])} = "monomodal"
        tform_type      (1,1) string {mustBeMember(tform_type, ["translation","rigid","similarity","affine"])} = "translation"
        max_step        (1,1) double {mustBePositive} = 1e-1
        min_step        (1,1) double {mustBePositive} = 1e-4
        max_gl_itn      (1,1) double {mustBePositive, mustBeInteger} = 50
        max_lo_itn      (1,:) double {mustBePositive, mustBeInteger} = 100
        iter_coeff      (1,1) double {mustBeInRange(iter_coeff, 0, 1)} = 0.5
        afs             (1,1) double {mustBeInRange(afs, 0.5, 3)} = 1.0
        vpl             (1,1) double {mustBePositive, mustBeInteger} = 3
        gl_interp       (1,1) string {mustBeMember(gl_interp, ["linear","cubic"])} = "linear"
        lo_interp       (1,1) string {mustBeMember(lo_interp, ["linear","cubic"])} = "cubic"

        % ================== one channel unique properties ================
        region_lbl      (1,1) string {mustBeMember(region_lbl, ["ORN","PN","LN","KC","MBON"])} = "ORN"
        gamma           (1,1) double {mustBeInRange(gamma, 0, 4)} = 2.0
        ocstdobj_th     (1,1) double {mustBeInRange(ocstdobj_th, 0, 5)} = 2.0
        ocscale_th      (1,1) double {mustBeNonnegative} = 3.0
        grid_unit       (1,1) string {mustBeMember(grid_unit, ["auto","1 1 1","2 2 1","3 3 1","4 4 1"])} = "auto"

        % ================ two channels unique properties ================
        autoctrst       (1,1) logical = false
        comacc          (1,1) double {mustBeMember(comacc, [512, 1024, 2048, 4096])} = 1024
        coregc          (1,1) double {mustBeInRange(coregc, 0, 100)} = 29
        gl_ds           (1,1) string {mustBeMember(gl_ds, ["auto","1X1","2X2","3X3"])} = "auto"
        strc_chl        (1,1) string {mustBeMember(strc_chl, ["r","g","b"])} = "r"
        func_chl        (1,1) string {mustBeMember(func_chl, ["r","g","b"])} = "g"

        % ==================== longterm properties ======================
        key_frames      (1,:) double {mustBePositive, mustBeInteger} = 1
        ltstdobj_th     (1,1) double {mustBeInRange(ltstdobj_th, 0, 5)} = 2.0
        ltscale_th      (1,1) double {mustBeNonnegative} = 3.0
        ltpct_ds        (1,1) string {mustBeMember(ltpct_ds, ["GA","RAND","NU"])} = "GA"
        ltvox_ds        (1,1) string {mustBeMember(ltvox_ds, ["1X1","2X2","3X3"])} = "2X2"
        ltpct_ds_arg    (1,1) double {mustBePositive} = 2.0
        ltpct_tol       (1,1) double {mustBePositive} = 1e-5
        lt_olr          (1,1) double {mustBeInRange(lt_olr, 0, 1)} = 0.1
        lt_iter_coeff   (1,1) double {mustBeInRange(lt_iter_coeff, 0, 1)} = 0.5
        max_ltpct_itn   (1,1) double {mustBePositive, mustBeInteger} = 50
        max_ltvox_itn   (1,:) double {mustBePositive, mustBeInteger} = 100
        max_lt_step     (1,1) double {mustBePositive} = 1e-2
        min_lt_step     (1,1) double {mustBePositive} = 1e-6
        lt_autoctrst    (1,1) logical = true
        lt_comacc       (1,1) double {mustBeMember(lt_comacc, [512, 1024, 2048, 4096])} = 1024
        lt_afs          (1,1) double {mustBeInRange(lt_afs, 0.5, 3)} = 1.0

        % =================== manual parameters ====================
        man_tform_type  (1,1) string {mustBeMember(man_tform_type, ...
                                    ["nonreflectivesimilarity","similarity","affine",...
                                    "projective","polynomial","pwl","lwm"])} ...
                                    = "nonreflectivesimilarity"
        degree          (1,1) double {mustBePositive, mustBeInteger} = 3
        nlwm            (1,1) double {mustBePositive, mustBeInteger} = 12
        man_interp      (1,1) string {mustBeMember(man_interp, ["linear","cubic"])} = "linear"
        edge_smooth     (1,1) logical = true      
    end

    properties(SetAccess=immutable, GetAccess=private, Hidden)
        reg_alg     (1,1) string {mustBeMember(reg_alg, ["OCREG","TCREG","MANREG","LTREG"])} = "TCREG"
        reg_mode    (1,1) string {mustBeMember(reg_mode, ["global","local","none"])} = "global"
        hardware    (1,1) string {mustBeMember(hardware, ["cpu", "cpu|gpu"])} = "cpu"
    end

    % properties(GetAccess={?RegisterWorker,?Register}, SetAccess = private, Dependent)
    %     Options
    % end

    properties(GetAccess=public, Dependent)
        Options
        Algorithm
        Mode
    end

    methods
        function this = regopt(regalg_, regmode_)
            %REGOPT A Constructor
            arguments
                regalg_   (1,1)   string {mustBeMember(regalg_, ...
                    ["OCREG","TCREG","MANREG","LTREG"])} = "TCREG"
                regmode_  (1,1)   string {mustBeMember(regmode_, ...
                    ["global","local","none"])} = "global"
            end

            this.reg_alg = regalg_;
            this.reg_mode = regmode_;
            this.hardware = get_hardware();
        end

        function r = get.Options(this)
            r = this.gen_option_();
        end

        function r = get.Algorithm(this)
            r = this.reg_alg;
        end

        function r = get.Mode(this)
            r = this.reg_mode;
        end
    end

    methods(Access=public, Hidden)
        function this = set(this, varargin)
            p = inputParser;
            % =========== not overloading parameters =============
            addParameter(p, 'RegModal',     this.reg_modal);
            addParameter(p, 'VPL',          this.vpl);
            addParameter(p, 'RL',           this.region_lbl);
            addParameter(p, 'Gamma',        this.gamma);
            addParameter(p, 'GridUnit',     this.grid_unit);
            addParameter(p, 'CoRegC',       this.coregc);
            addParameter(p, 'DS',           this.gl_ds);
            addParameter(p, 'SC',           this.strc_chl);
            addParameter(p, 'FC',           this.func_chl);
            addParameter(p, 'KeyFrames',    this.key_frames);
            addParameter(p, 'PCTDS',        this.ltpct_ds);
            addParameter(p, 'PCTDSArg',     this.ltpct_ds_arg);
            addParameter(p, 'VoxDS',        this.ltvox_ds);
            addParameter(p, 'TolPCT',       this.ltpct_tol);
            addParameter(p, 'Outlier',      this.lt_olr);
            addParameter(p, 'MaxPCTIterN',  this.max_ltpct_itn);
            addParameter(p, 'MaxVoxIterN',  this.max_ltvox_itn);
            addParameter(p, 'Degree',       this.degree);
            addParameter(p, 'Nlwm',         this.nlwm);
            addParameter(p, 'EdgeSmooth',   this.edge_smooth);

            % =============== overloading parameters =============
            switch this.reg_alg
                case "OCREG"
                    switch this.reg_mode
                        case "global"
                            addParameter(p, 'TformType',    this.tform_type);
                            addParameter(p, 'MaxStep',      this.max_step);
                            addParameter(p, 'MinStep',      this.min_step);
                            addParameter(p, 'MaxIterN',     this.max_gl_itn);
                            addParameter(p, 'IterCoeff',    this.iter_coeff);
                            addParameter(p, 'Interp',       this.gl_interp);
                            addParameter(p, 'ThFG',         this.ocstdobj_th);
                            addParameter(p, 'ThScale',      this.ocscale_th);
                        case "local"
                            addParameter(p, 'MaxIterN',     this.max_lo_itn);
                            addParameter(p, 'AFS',          this.afs);
                            addParameter(p, 'Interp',       this.lo_interp);
                    end
                case "TCREG"
                    switch this.reg_mode
                        case "global"
                            addParameter(p, 'TformType',    this.tform_type);
                            addParameter(p, 'MaxStep',      this.max_step);
                            addParameter(p, 'MinStep',      this.min_step);
                            addParameter(p, 'MaxIterN',     this.max_gl_itn);
                            addParameter(p, 'IterCoeff',    this.iter_coeff);
                            addParameter(p, 'Interp',       this.gl_interp);
                        case "local"
                            addParameter(p, 'MaxIterN',     this.max_lo_itn);
                            addParameter(p, 'AFS',          this.afs);
                            addParameter(p, 'Interp',       this.lo_interp);
                            addParameter(p, 'AutoContrast', this.autoctrst);
                            addParameter(p, 'ComAcc',       this.comacc);
                    end
                case "MANREG"
                    addParameter(p, 'TformType',    this.man_tform_type);
                    addParameter(p, 'Interp',       this.man_interp);
                case "LTREG"
                    addParameter(p, 'ThFG',         this.ltstdobj_th);
                    addParameter(p, 'ThScale',      this.ltscale_th);
                    addParameter(p, 'MaxStep',      this.max_lt_step);
                    addParameter(p, 'MinStep',      this.min_lt_step);
                    addParameter(p, 'IterCoeff',    this.lt_iter_coeff);
                    addParameter(p, 'AFS',          this.lt_afs);
                    addParameter(p, 'AutoContrast', this.lt_autoctrst);
                    addParameter(p, 'ComAcc',       this.lt_comacc);
                otherwise
            end

            parse(p, varargin{:});

            this = set_option_(this, p.Results);
        end
    end


    methods(Access=private)
        function r = gen_option_(this)
            switch this.reg_alg
                case "OCREG"
                    switch this.reg_mode
                        case "global"
                            r = struct("RegModal",  this.reg_modal, ...
                                       "TformType", this.tform_type, ...
                                       "MaxStep",   this.max_step, ...
                                       "MinStep",   this.min_step, ...
                                       "MaxIterN",  this.max_gl_itn, ...
                                       "IterCoeff", this.iter_coeff, ...
                                       "VPL",       this.vpl, ...
                                       "Interp",    this.gl_interp, ...
                                       "RL",        this.region_lbl, ...
                                       "Gamma",     this.gamma, ...
                                       "ThFG",      this.ocstdobj_th, ...
                                       "ThScale",   this.ocscale_th, ...
                                       "Hardware",  this.hardware);
                        case "local"
                            r = struct("MaxIterN",  this.max_lo_itn, ...
                                       "AFS",       this.afs, ...
                                       "VPL",       this.vpl, ...
                                       "Interp",    this.lo_interp, ...
                                       "RL",        this.region_lbl, ...
                                       "Gamma",     this.gamma, ...
                                       "GridUnit",  this.grid_unit, ...
                                       "Hardware",  this.hardware);
                        otherwise
                    end
                case "TCREG"
                    switch this.reg_mode
                        case "global"
                            r = struct("RegModal",      this.reg_modal, ...
                                       "TformType",     this.tform_type, ...
                                       "MaxStep",       this.max_step, ...
                                       "MinStep",       this.min_step, ...
                                       "MaxIterN",      this.max_gl_itn, ...
                                       "IterCoeff",     this.iter_coeff, ...
                                       "VPL",           this.vpl, ...
                                       "Interp",        this.gl_interp, ...
                                       "CoRegC",        this.coregc, ...
                                       "DS",            this.gl_ds, ...
                                       "SC",            this.strc_chl, ...
                                       "FC",            this.func_chl, ...
                                       "Hardware",      this.hardware);
                        case "local"
                            r = struct("MaxIterN",      this.max_lo_itn, ...
                                       "AFS",           this.afs, ...
                                       "VPL",           this.vpl, ...
                                       "Interp",        this.lo_interp, ...
                                       "AutoContrast",  this.autoctrst, ...
                                       "ComAcc",        this.comacc, ...
                                       "Hardware",      this.hardware);
                        otherwise
                    end
                case "MANREG"
                    r = struct("TformType",     this.man_tform_type, ...
                               "Degree",        this.degree, ...
                               "Nlwm",          this.nlwm, ...
                               "Interp",        this.man_interp, ...
                               "EdgeSmooth",    this.edge_smooth, ...
                               "Hardware",      this.hardware);
                case "LTREG"
                    r = struct("KeyFrames",     this.key_frames, ...
                               "ThFG",          this.ltstdobj_th, ...
                               "ThScale",       this.ltscale_th, ...
                               "PCTDS",         this.ltpct_ds, ...
                               "PCTDSArg",      this.ltpct_ds_arg, ...
                               "VoxDS",         this.ltvox_ds, ...
                               "TolPCT",        this.ltpct_tol, ...
                               "Outlier",       this.lt_olr, ...
                               "IterCoeff",     this.lt_iter_coeff, ...
                               "MaxPCTIterN",   this.max_ltpct_itn, ...
                               "MaxVoxIterN",   this.max_ltvox_itn, ...
                               "MaxStep",       this.max_lt_step, ...
                               "MinStep",       this.min_lt_step, ...
                               "AutoContrast",  this.lt_autoctrst, ...
                               "ComAcc",        this.lt_comacc, ...
                               "AFS",           this.lt_afs, ...
                               "Hardware",      this.hardware);
                otherwise

            end
        end

        function this = set_option_(this, r_)
            % r_ is parser results struct, inverse the mapping
            switch this.reg_alg
                case "OCREG"
                    switch this.reg_mode
                        case "global"
                            this.reg_modal = r_.RegModal;
                            this.tform_type = r_.TformType;
                            this.max_step = r_.MaxStep;
                            this.min_step = r_.MinStep;
                            this.max_gl_itn = r_.MaxIterN;
                            this.iter_coeff = r_.IterCoeff;
                            this.vpl = r_.VPL;
                            this.gl_interp = r_.Interp;
                            this.region_lbl = r_.RL;
                            this.gamma = r_.Gamma;
                            this.ocstdobj_th = r_.ThFG;
                            this.ocscale_th = r_.ThScale;
                        case "local"
                            this.max_lo_itn = r_.MaxIterN;
                            this.afs = r_.AFS;
                            this.vpl = r_.VPL;
                            this.lo_interp = r_.Interp;
                            this.region_lbl = r_.RL;
                            this.gamma = r_.Gamma;
                            this.grid_unit = r_.GridUnit;
                        otherwise
                    end
                case "TCREG"
                    switch this.reg_mode
                        case "global"
                            this.reg_modal = r_.RegModal;
                            this.tform_type = r_.TformType;
                            this.max_step = r_.MaxStep;
                            this.min_step = r_.MinStep;
                            this.max_gl_itn = r_.MaxIterN;
                            this.iter_coeff = r_.IterCoeff;
                            this.vpl = r_.VPL;
                            this.gl_interp = r_.Interp;
                            this.coregc = r_.CoRegC;
                            this.gl_ds = r_.DS;
                            this.strc_chl = r_.SC;
                            this.func_chl = r_.FS;
                        case "local"
                            this.max_lo_itn = r_.MaxIterN;
                            this.afs = r_.AFS;
                            this.vpl = r_.VPL;
                            this.lo_interp = r_.Interp;
                            this.autoctrst = r_.AutoContrast;
                            this.comacc = r_.ComAcc;
                    end
                case "MANREG"
                    this.man_tform_type = r_.TformType;
                    this.degree = r_.Degree;
                    this.nlwm = r_.Nlwm;
                    this.man_interp = r_.Interp;
                    this.edge_smooth = r_.EdgeSmooth;
                case "LTREG"
                    this.key_frames = r_.KeyFrames;
                    this.ltstdobj_th = r_.ThFG;
                    this.ltscale_th = r_.ThScale;
                    this.ltpct_ds = r_.PCTDS;
                    this.ltpct_ds_arg = r_.pctDSArg;
                    this.ltvox_ds = r_.VoxDS;
                    this.ltpct_tol = r_.TolPCT;
                    this.lt_olr = r_.Outlier;
                    this.lt_iter_coeff = r_IterCoeff;
                    this.max_ltpct_itn = r_.MaxPCTIterN;
                    this.max_ltvox_itn = r_.MaxVoxIterN;
                    this.max_lt_step = r_.MaxStep;
                    this.min_lt_step = r_.MinStep;
                    this.lt_autoctrst = r_.AutoContrast;
                    this.lt_comacc = r_.ComAcc;
                    this.lt_afs = r_.AFS;
                otherwise
            end
        end
    end
end

function r = get_hardware()
if gpuDeviceCount > 0
    r = "cpu|gpu";
else
    r = "cpu";
end
end