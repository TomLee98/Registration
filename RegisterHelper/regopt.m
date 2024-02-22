classdef regopt
    %REGOPT This class is registration option defination, which is packaged
    %value class
    
    properties(Access=private, Hidden)
        % ================ all options common properties =================
        reg_alg         (1,1) string {mustBeMember(reg_alg, ["OCREG","TCREG","MANREG","LTREG"])} = "TCREG"
        reg_mode        (1,1) string {mustBeMember(reg_mode, ["global","local"])} = "global"
        sub_reg_alg     (1,1) string {mustBeMember(sub_reg_alg, ["usual","advanced"])} = "usual"

        % ============= one/two channel(s) common properties ==============
        reg_modal       (1,1) string {mustBeMember(reg_modal, ["multimodal", "monomodal"])} = "monomodal"
        tform_type      (1,1) string {mustBeMember(tform_type, ["translation","rigid","affine"])} = "translation"
        step_max        (1,1) double {mustBePositive} = 1e-1
        step_min        (1,1) double {mustBePositive} = 1e-4
        gl_itn_max      (1,1) double {mustBePositive, mustBeInteger} = 50
        lo_itn_max      (1,:) double {mustBePositive, mustBeInteger} = 100
        iter_coeff      (1,1) double {mustBeInRange(iter_coeff, 0, 1)} = 0.5
        afs             (1,1) double {mustBeInRange(afs, 0.5, 3)} = 1.0
        vpl             (1,1) double {mustBePositive, mustBeInteger} = 3
        gl_interp       (1,1) string {mustBeMember(gl_interp, ["linear","cubic"])} = "linear"
        lo_interp       (1,1) string {mustBeMember(lo_interp, ["linear","cubic"])} = "linear"
        mfilter         (1,3) double {mustBeNonnegative, mustBeInteger} = [3,3,3]
        open_operator   (1,3) double {mustBeNonnegative, mustBeInteger} = [5,5,2]
        gfilter         (1,3) double {mustBeNonnegative, mustBeInteger} = [3,3,3]
        gamma           (1,1) double {mustBeInRange(gamma, 0, 4)} = 1.0
        zopt_shift_max  (1,1) double {mustBeNonnegative} = 2
        zopt_tol        (1,1) double {mustBeInRange(zopt_tol, 0, 1)} = 1e-3

        % ================== one channel unique properties ================
        region_lbl      (1,1) string {mustBeMember(region_lbl, ["ORN","PN","LN","KC","MBON"])} = "ORN"
        ocstdobj_th     (1,1) double {mustBeInRange(ocstdobj_th, 0, 5)} = 2.0
        ocscale_th      (1,1) double {mustBeNonnegative} = 3.0
        grid_unit       (1,1) string {mustBeMember(grid_unit, ["auto","1 1 1","2 2 1","3 3 1","4 4 1"])} = "auto"

        % ================ two channels unique properties ================
        img_rehist      (1,1) logical = false
        repacc          (1,1) double {mustBeMember(repacc, [2048, 4096, 8192])} = 2048
        grid_regulation (1,1) double {mustBeInRange(grid_regulation, 0, 5)} = 0.11
        grid_spacing    (1,3) double {mustBePositive} = [4,4,4]
        ds              (1,1) string {mustBeMember(ds, ["auto","1X1","2X2","3X3"])} = "auto"
        strc_chl        (1,1) string {mustBeMember(strc_chl, ["r","g","b"])} = "r"
        func_chl        (1,1) string {mustBeMember(func_chl, ["r","g","b"])} = "g"

        % ==================== longterm properties ======================
        lt_keyframes    (1,:) double {mustBePositive, mustBeInteger} = 1
        lt_autokey      (1,1) logical = true
        lt_stdobj_th    (1,1) double {mustBeInRange(lt_stdobj_th, 0, 5)} = 2.0
        lt_scale_th     (1,1) double {mustBeNonnegative} = 3.0
        lt_ds           (1,1) string {mustBeMember(lt_ds, ["GA","RAND","NU"])} = "GA"
        lt_ds_arg       (1,1) double {mustBePositive} = 2.0
        lt_tol          (1,1) double {mustBePositive} = 1e-5
        lt_olr          (1,1) double {mustBeInRange(lt_olr, 0, 1)} = 0.1
        lt_iter_coeff   (1,1) double {mustBeInRange(lt_iter_coeff, 0, 1)} = 0.5
        lt_itn_max      (1,1) double {mustBePositive, mustBeInteger} = 50
        lt_step_max     (1,1) double {mustBePositive} = 1e-2
        lt_step_min     (1,1) double {mustBePositive} = 1e-6

        % =================== manual parameters ====================
        m_tform_type     (1,1) string {mustBeMember(m_tform_type, ...
                                    ["translation","rigid","poly","pwl"])} ...
                                    = "translation"
        m_degree        (1,1) double {mustBeMember(m_degree, [2,3,4])} = 2
        m_dview         (1,1) string {mustBeMember(m_dview, ["XY","ZX","YZ"])} = "XY"
        m_proj          (1,1) string {mustBeMember(m_proj, ["max","min","median","mean"])} = "max"
        m_interp        (1,1) string {mustBeMember(m_interp, ["linear","cubic"])} = "linear"
        m_rs            (1,3) double {mustBePositive} = [1,1,1]
        m_isometric     (1,1) logical = false
    end

    properties(SetAccess=immutable, GetAccess=private, Hidden)
        hardware    (1,1) string {mustBeMember(hardware, ["cpu", "cpu|gpu"])} = "cpu"
    end

    properties(GetAccess=public, Dependent)
        Algorithm
        Mode
        Options
        SubAlgorithm
    end

    methods
        function this = regopt(regalg_, regmode_)
            %REGOPT A Constructor
            arguments
                regalg_   (1,1)   string {mustBeMember(regalg_, ...
                    ["OCREG","TCREG","MANREG","LTREG"])} = "TCREG"
                regmode_  (1,1)   string {mustBeMember(regmode_, ...
                    ["global","local"])} = "global"
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

        function this = set.Algorithm(this, r_)
            arguments
                this
                r_   (1,1)   string {mustBeMember(r_, ...
                    ["OCREG","TCREG","MANREG","LTREG"])} = "TCREG"
            end

            this.reg_alg = r_;
        end

        function r = get.Mode(this)
            r = this.reg_mode;
        end

        function this = set.Mode(this, r_)
            arguments
                this
                r_  (1,1)   string {mustBeMember(r_, ["global","local"])} ...
                    = "global"
            end

            this.reg_mode = r_;
        end

        function r = get.SubAlgorithm(this)
            r = this.sub_reg_alg;
        end
        
        function this = set.SubAlgorithm(this, r_)
            arguments
                this
                r_  (1,1)   string {mustBeMember(r_, ["usual","advanced"])} ...
                    = "usual"
            end

            this.sub_reg_alg = r_;
        end
    end

    methods(Access=public, Hidden)
        function this = set(this, varargin)
            p = inputParser;
            % =========== not overloading parameters =============
            addParameter(p, 'RegModal',         this.reg_modal);
            addParameter(p, 'MedianFilter',     this.mfilter);
            addParameter(p, 'OpenOperator',     this.open_operator);
            addParameter(p, 'GaussianFilter',   this.gfilter);
            addParameter(p, 'MaxZOptShift',     this.zopt_shift_max);
            addParameter(p, 'TolZOpt',          this.zopt_tol);
            addParameter(p, 'VPL',              this.vpl);
            addParameter(p, 'RL',               this.region_lbl);
            addParameter(p, 'Gamma',            this.gamma);
            addParameter(p, 'GridUnit',         this.grid_unit);
            addParameter(p, 'SC',               this.strc_chl);
            addParameter(p, 'FC',               this.func_chl);
            addParameter(p, 'KeyFrames',        this.lt_keyframes);
            addParameter(p, 'AutoKeyFrame',     this.lt_autokey);
            addParameter(p, 'Tol',              this.lt_tol);
            addParameter(p, 'Outlier',          this.lt_olr);

            % =============== overloading parameters =============
            switch this.reg_alg
                case "OCREG"
                    switch this.reg_mode
                        case "global"
                            addParameter(p, 'TformType',    this.tform_type);
                            addParameter(p, 'MaxStep',      this.step_max);
                            addParameter(p, 'MinStep',      this.step_min);
                            addParameter(p, 'MaxIterN',     this.gl_itn_max);
                            addParameter(p, 'IterCoeff',    this.iter_coeff);
                            addParameter(p, 'Interp',       this.gl_interp);
                            addParameter(p, 'ThFG',         this.ocstdobj_th);
                            addParameter(p, 'ThScale',      this.ocscale_th);
                            addParameter(p, 'DS',           this.ds);
                        case "local"
                            addParameter(p, 'MaxIterN',     this.lo_itn_max);
                            addParameter(p, 'AFS',          this.afs);
                            addParameter(p, 'Interp',       this.lo_interp);
                        otherwise
                    end
                case "TCREG"
                    switch this.reg_mode
                        case "global"
                            addParameter(p, 'TformType',    this.tform_type);
                            addParameter(p, 'MaxStep',      this.step_max);
                            addParameter(p, 'MinStep',      this.step_min);
                            addParameter(p, 'MaxIterN',     this.gl_itn_max);
                            addParameter(p, 'IterCoeff',    this.iter_coeff);
                            addParameter(p, 'Interp',       this.gl_interp);
                            addParameter(p, 'DS',           this.ds);
                        case "local"
                            addParameter(p, 'MaxIterN',     this.lo_itn_max);
                            addParameter(p, 'AFS',          this.afs);
                            addParameter(p, 'Interp',       this.lo_interp);
                            addParameter(p, 'ImageRehist',  this.img_rehist);
                            addParameter(p, 'RepAcc',       this.repacc);
                            addParameter(p, 'GR',           this.grid_regulation);
                            addParameter(p, 'GS',           this.grid_spacing);
                        otherwise
                    end
                case "MANREG"
                    addParameter(p, 'TformType',    this.m_tform_type);
                    addParameter(p, 'Interp',       this.m_interp);
                    addParameter(p, 'Degree',       this.m_degree);
                    addParameter(p, 'DView',        this.m_dview);
                    addParameter(p, 'Projection',   this.m_proj);
                    addParameter(p, 'Resampling',   this.m_rs);
                    addParameter(p, 'Isometric',    this.m_isometric);
                case "LTREG"
                    addParameter(p, 'TformType',    this.tform_type);
                    addParameter(p, 'MaxIterN',     this.lt_itn_max);
                    addParameter(p, 'ThFG',         this.lt_stdobj_th);
                    addParameter(p, 'ThScale',      this.lt_scale_th);
                    addParameter(p, 'MaxStep',      this.lt_step_max);
                    addParameter(p, 'MinStep',      this.lt_step_min);
                    addParameter(p, 'IterCoeff',    this.lt_iter_coeff);
                    addParameter(p, 'DS',           this.lt_ds);
                    addParameter(p, 'DSArg',        this.lt_ds_arg);
                otherwise
            end

            parse(p, varargin{:});

            this = this.set_option_(p.Results);
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
                                       "MaxStep",   this.step_max, ...
                                       "MinStep",   this.step_min, ...
                                       "MaxIterN",  this.gl_itn_max, ...
                                       "IterCoeff", this.iter_coeff, ...
                                       "VPL",       this.vpl, ...
                                       "Interp",    this.gl_interp, ...
                                       "RL",        this.region_lbl, ...
                                       "Gamma",     this.gamma, ...
                                       "ThFG",      this.ocstdobj_th, ...
                                       "ThScale",   this.ocscale_th, ...
                                       "Hardware",  this.hardware);
                        case "local"
                            r = struct("MaxIterN",  this.lo_itn_max, ...
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
                                       "MedianFilter",  this.mfilter, ...
                                       "OpenOperator",  this.open_operator, ...
                                       "GaussianFilter",this.gfilter, ...
                                       "MaxZOptShift",  this.zopt_shift_max, ...
                                       "TolZOpt",       this.zopt_tol, ...
                                       "Gamma",         this.gamma, ...
                                       "MaxStep",       this.step_max, ...
                                       "MinStep",       this.step_min, ...
                                       "MaxIterN",      this.gl_itn_max, ...
                                       "IterCoeff",     this.iter_coeff, ...
                                       "VPL",           this.vpl, ...
                                       "Interp",        this.gl_interp, ...
                                       "DS",            this.ds, ...
                                       "SC",            this.strc_chl, ...
                                       "FC",            this.func_chl, ...
                                       "Hardware",      this.hardware);
                        case "local"
                            r = struct("MaxIterN",      this.lo_itn_max, ...
                                       "AFS",           this.afs, ...
                                       "GR",            this.grid_regulation, ...
                                       "GS",            this.grid_spacing, ...
                                       "VPL",           this.vpl, ...
                                       "Interp",        this.lo_interp, ...
                                       "ImageRehist",   this.img_rehist, ...
                                       "RepAcc",        this.repacc, ...
                                       "SC",            this.strc_chl, ...
                                       "FC",            this.func_chl, ...
                                       "Hardware",      this.hardware);
                        otherwise
                    end
                case "MANREG"
                    r = struct("TformType",     this.m_tform_type, ...
                               "Degree",        this.m_degree, ...
                               "Interp",        this.m_interp, ...
                               "DView",         this.m_dview, ...
                               "Projection",    this.m_proj, ...
                               "Resampling",    this.m_rs, ...
                               "Isometric",     this.m_isometric, ...
                               "SC",            this.strc_chl, ...
                               "FC",            this.func_chl, ...
                               "Hardware",      this.hardware);
                case "LTREG"
                    r = struct("TformType",     this.tform_type, ...
                               "KeyFrames",     this.lt_keyframes, ...
                               "AutoKeyFrame",  this.lt_autokey, ...
                               "ThFG",          this.lt_stdobj_th, ...
                               "ThScale",       this.lt_scale_th, ...
                               "DS",            this.lt_ds, ...
                               "DSArg",         this.lt_ds_arg, ...
                               "Tol",           this.lt_tol, ...
                               "Outlier",       this.lt_olr, ...
                               "IterCoeff",     this.lt_iter_coeff, ...
                               "MaxIterN",      this.lt_itn_max, ...
                               "MaxStep",       this.lt_step_max, ...
                               "MinStep",       this.lt_step_min, ...
                               "SC",            this.strc_chl, ...
                               "FC",            this.func_chl, ...
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
                            this.step_max = r_.MaxStep;
                            this.step_min = r_.MinStep;
                            this.gl_itn_max = r_.MaxIterN;
                            this.iter_coeff = r_.IterCoeff;
                            this.vpl = r_.VPL;
                            this.gl_interp = r_.Interp;
                            this.region_lbl = r_.RL;
                            this.gamma = r_.Gamma;
                            this.ocstdobj_th = r_.ThFG;
                            this.ocscale_th = r_.ThScale;
                        case "local"
                            this.lo_itn_max = r_.MaxIterN;
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
                            this.mfilter = r_.MedianFilter;
                            this.open_operator = r_.OpenOperator;
                            this.gfilter = r_.GaussianFilter;
                            this.zopt_shift_max = r_.MaxZOptShift;
                            this.zopt_tol = r_.TolZOpt;
                            this.gamma = r_.Gamma;
                            this.step_max = r_.MaxStep;
                            this.step_min = r_.MinStep;
                            this.gl_itn_max = r_.MaxIterN;
                            this.iter_coeff = r_.IterCoeff;
                            this.vpl = r_.VPL;
                            this.gl_interp = r_.Interp;
                            this.ds = r_.DS;
                            this.strc_chl = r_.SC;
                            this.func_chl = r_.FC;
                        case "local"
                            this.lo_itn_max = r_.MaxIterN;
                            this.afs = r_.AFS;
                            this.grid_regulation = r_.GR;
                            this.grid_spacing = r_.GS;
                            this.vpl = r_.VPL;
                            this.lo_interp = r_.Interp;
                            this.img_rehist = r_.ImageRehist;
                            this.repacc = r_.RepAcc;
                            this.strc_chl = r_.SC;
                            this.func_chl = r_.FC;
                    end
                case "MANREG"
                    this.m_tform_type = r_.TformType;
                    this.m_degree = r_.Degree;
                    this.m_interp = r_.Interp;
                    this.m_dview = r_.DView;
                    this.m_proj = r_.Projection;
                    this.m_rs = r_.Resampling;
                    this.m_isometric = r_.Isometric;
                    this.strc_chl = r_.SC;
                    this.func_chl = r_.FC;
                case "LTREG"
                    this.tform_type = r_.TformType;
                    this.lt_keyframes = r_.KeyFrames;
                    this.lt_autokey = r_.AutoKeyFrame;
                    this.lt_stdobj_th = r_.ThFG;
                    this.lt_scale_th = r_.ThScale;
                    this.lt_ds = r_.DS;
                    this.lt_ds_arg = r_.DSArg;
                    this.lt_tol = r_.Tol;
                    this.lt_olr = r_.Outlier;
                    this.lt_iter_coeff = r_.IterCoeff;
                    this.lt_itn_max = r_.MaxIterN;
                    this.lt_step_max = r_.MaxStep;
                    this.lt_step_min = r_.MinStep;
                    this.strc_chl = r_.SC;
                    this.func_chl = r_.FC;
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