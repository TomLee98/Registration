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
        step_max        (1,1) double {mustBePositive} = 1e-1
        step_min        (1,1) double {mustBePositive} = 1e-4
        gl_itn_max      (1,1) double {mustBePositive, mustBeInteger} = 50
        gl_interp       (1,1) string {mustBeMember(gl_interp, ["linear","cubic"])} = "linear"
        iter_coeff      (1,1) double {mustBePositive} = 0.5
        gl_vpl          (1,1) double {mustBePositive, mustBeInteger} = 3
        zopt_shift_max  (1,1) double {mustBeNonnegative} = 2
        zopt_tol        (1,1) double {mustBeInRange(zopt_tol, 0, 1)} = 1e-3
        strc_chl        (1,1) string {mustBeMember(strc_chl, ["r","g","b",""])} = ""
        func_chl        (1,1) string {mustBeMember(func_chl, ["r","g","b",""])} = ""

        % ================== one channel unique properties ================
        % [global]:
        oc_tform_type   (1,1) string {mustBeMember(oc_tform_type, ["translation","rigid","affine"])} = "translation"
        oc_coarse_alg   (1,1) string {mustBeMember(oc_coarse_alg, ["mmt","pcorr","none"])} = "mmt"
        oc_coarse_args  (1,1) struct = struct("Filter", "median", "VT", 1000, "Radius", 3)
        oc_mfilter      (1,3) double {mustBeNonnegative, mustBeInteger} = [7,7,3]
        oc_gfilter      (1,3) double {mustBeNonnegative, mustBeInteger} = [1,1,1]

        % ================ two channels unique properties ================
        % [global]:
        tc_tform_type   (1,1) string {mustBeMember(tc_tform_type, ["translation","rigid","affine"])} = "translation"
        tc_coarse_alg   (1,1) string {mustBeMember(tc_coarse_alg, ["mmt","pcorr","fpp","none"])} = "mmt"
        tc_coarse_args  (1,1) struct = struct("Operator", "SIFT", "QT", 0.0133, "NumOctave", 3)
        tc_mfilter      (1,3) double {mustBeNonnegative, mustBeInteger} = [3,3,3]
        tc_gfilter      (1,3) double {mustBeNonnegative, mustBeInteger} = [3,3,3]
        dfilter         (1,4) double {mustBeInteger} = [3,115,65535,1000]
        dfilter_enh     (1,1) logical = false
        area_mask       (1,1) logical = false
        gamma           (1,1) double {mustBeInRange(gamma, 0, 4)} = 1.0
        ds              (1,1) string {mustBeMember(ds, ["auto","1X1","2X2","3X3"])} = "auto"
        % [local]:
        lo_itn_max      (1,:) double {mustBePositive, mustBeInteger} = 100
        lo_afs          (1,1) double {mustBeInRange(lo_afs, 0.5, 3)} = 1.0
        lo_interp       (1,1) string {mustBeMember(lo_interp, ["linear","cubic"])} = "linear"
        img_rehist      (1,1) logical = false
        dm_vpl          (1,1) double {mustBePositive, mustBeInteger} = 3    % vpl for imregdemons
        df_vpl          (1,1) double {mustBePositive, mustBeInteger} = 3    % vpl for imregdeform
        repacc          (1,1) double {mustBeMember(repacc, [2048, 4096, 8192])} = 2048
        grid_regulation (1,1) double {mustBeInRange(grid_regulation, 0, 5)} = 0.11
        grid_spacing    (1,3) double {mustBePositive} = [4,4,4]

        % ==================== longterm properties ======================
        lt_keyframes    (:,1) double {mustBePositive, mustBeInteger} = 1
        lt_autokey      (1,1) logical = true
        lt_autotpl      (1,1) logical  = true
        lt_tgridminmax  (1,2) double {mustBePositive, mustBeInteger} = [10, 40] % typical 60*[10,40]=[600, 2400] volumes is best range
        lt_regchain     (1,1) regchain = regchain.empty()
        lt_iter_coeff   (1,1) double {mustBeInRange(lt_iter_coeff, 0, 1)} = 0.5
        lt_itn_max      (1,1) double {mustBePositive, mustBeInteger} = 100
        lt_step_max     (1,1) double {mustBePositive} = 1e-2
        lt_step_min     (1,1) double {mustBePositive} = 1e-6
        lt_interp       (1,1) string {mustBeMember(lt_interp, ["linear","cubic"])} = "linear"
        
        % =================== manual parameters ====================
        m_tform_type     (1,1) string {mustBeMember(m_tform_type, ...
                                    ["translation","rigid","affine","poly","pwl"])} ...
                                    = "translation"
        m_degree        (1,1) double {mustBeMember(m_degree, [2,3,4])} = 2
        m_dview         (1,1) string {mustBeMember(m_dview, ["XY","ZX","YZ"])} = "XY"
        m_proj          (1,1) string {mustBeMember(m_proj, ["max","min","median","mean"])} = "max"
        m_interp        (1,1) string {mustBeMember(m_interp, ["linear","cubic"])} = "linear"
        m_rs            (1,3) double {mustBePositive} = [1,1,1]
        m_isometric     (1,1) logical = false
        m_tfmat         (4,4) double {mustBeFinite} = eye(4)
        m_tfmat_enable  (1,1) logical = false
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
            p.StructExpand = false;     % allow struct as parameter

            % =========== not overloading parameters =============
            addParameter(p, 'RegModal',     this.reg_modal);
            addParameter(p, 'MaxZOptShift', this.zopt_shift_max);
            addParameter(p, 'TolZOpt',      this.zopt_tol);
            addParameter(p, 'glVPL',        this.gl_vpl);
            addParameter(p, 'DS',           this.ds);
            addParameter(p, 'SC',           this.strc_chl);
            addParameter(p, 'FC',           this.func_chl);

            % =========== overloading parameters (shared name) ===========
            switch this.reg_alg
                case "OCREG"
                    switch this.reg_mode
                        case "global"
                            addParameter(p, 'TformType',    this.oc_tform_type);
                            addParameter(p, 'MedianFilter',     this.oc_mfilter);
                            addParameter(p, 'GaussianFilter',   this.oc_gfilter);
                            addParameter(p, 'MaxStep',      this.step_max);
                            addParameter(p, 'MinStep',      this.step_min);
                            addParameter(p, 'MaxIterN',     this.gl_itn_max);
                            addParameter(p, 'IterCoeff',    this.iter_coeff);
                            addParameter(p, 'Interp',       this.gl_interp);
                            addParameter(p, 'CoarseAlg',    this.oc_coarse_alg);
                            addParameter(p, 'CoarseArgs',   this.oc_coarse_args);
                        otherwise
                    end
                case "TCREG"
                    switch this.reg_mode
                        case "global"
                            addParameter(p, 'Gamma',        this.gamma);
                            addParameter(p, 'MedianFilter',     this.tc_mfilter);
                            addParameter(p, 'GaussianFilter',   this.tc_gfilter);
                            addParameter(p, 'DilateFilter',     this.dfilter);
                            addParameter(p, 'DilateFilterEnh',  this.dfilter_enh);
                            addParameter(p, 'AreaMask',     this.area_mask);
                            addParameter(p, 'TformType',    this.tc_tform_type);
                            addParameter(p, 'MaxStep',      this.step_max);
                            addParameter(p, 'MinStep',      this.step_min);
                            addParameter(p, 'MaxIterN',     this.gl_itn_max);
                            addParameter(p, 'IterCoeff',    this.iter_coeff);
                            addParameter(p, 'Interp',       this.gl_interp);
                            addParameter(p, 'CoarseAlg',    this.tc_coarse_alg);
                            addParameter(p, 'CoarseArgs',   this.tc_coarse_args);
                        case "local"
                            addParameter(p, 'MaxIterN',     this.lo_itn_max);
                            addParameter(p, 'AFS',          this.lo_afs);
                            addParameter(p, 'Interp',       this.lo_interp);
                            addParameter(p, 'ImageRehist',  this.img_rehist);
                            addParameter(p, 'RepAcc',       this.repacc);
                            addParameter(p, 'GR',           this.grid_regulation);
                            addParameter(p, 'GS',           this.grid_spacing);
                            addParameter(p, 'dmVPL',        this.dm_vpl);
                            addParameter(p, 'dfVPL',        this.df_vpl);
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
                    addParameter(p, 'TFMatrix',     this.m_tfmat);
                    addParameter(p, 'TFEnable',     this.m_tfmat_enable);
                case "LTREG"
                    addParameter(p, 'MaxIterN',     this.lt_itn_max);
                    addParameter(p, 'MaxStep',      this.lt_step_max);
                    addParameter(p, 'MinStep',      this.lt_step_min);
                    addParameter(p, 'IterCoeff',    this.lt_iter_coeff);
                    addParameter(p, 'Interp',       this.lt_interp);
                    addParameter(p, 'Keyframes',    this.lt_keyframes);
                    addParameter(p, 'AutoKeyframe', this.lt_autokey);
                    addParameter(p, 'AutoTemplate', this.lt_autotpl);
                    addParameter(p, 'TGridMinMax',  this.lt_tgridminmax);
                    addParameter(p, 'RegChain',     this.lt_regchain);
                    addParameter(p, 'DilateFilter', this.dfilter);
                otherwise
            end

            parse(p, varargin{:});

            this = this.set_option_(p.Results);
        end
    end


    methods(Access=private)
        function r = gen_option_(this)
            % this is a filter function
            switch this.reg_alg
                case "OCREG"
                    switch this.reg_mode
                        case "global"
                            r = struct("RegModal",      this.reg_modal, ...
                                       "TformType",     this.oc_tform_type, ...
                                       "GaussianFilter",this.oc_gfilter, ...
                                       "MedianFilter",  this.oc_mfilter, ...
                                       "MaxStep",       this.step_max, ...
                                       "MinStep",       this.step_min, ...
                                       "MaxIterN",      this.gl_itn_max, ...
                                       "CoarseAlg",     this.oc_coarse_alg, ...
                                       "CoarseArgs",    this.oc_coarse_args, ...
                                       "MaxZOptShift",  this.zopt_shift_max, ...
                                       "TolZOpt",       this.zopt_tol, ...
                                       "IterCoeff",     this.iter_coeff, ...
                                       "glVPL",         this.gl_vpl, ...
                                       "Interp",        this.gl_interp, ...
                                       "DS",            this.ds, ...
                                       "SC",            this.strc_chl, ...
                                       "FC",            this.func_chl, ...
                                       "Hardware",      this.hardware);
                        otherwise
                    end
                case "TCREG"
                    switch this.reg_mode
                        case "global"
                            r = struct("RegModal",      this.reg_modal, ...
                                       "AreaMask",      this.area_mask, ...
                                       "TformType",     this.tc_tform_type, ...
                                       "MedianFilter",  this.tc_mfilter, ...
                                       "GaussianFilter",this.tc_gfilter, ...
                                       "DilateFilter",  this.dfilter, ...
                                       "DilateFilterEnh",this.dfilter_enh, ...
                                       "MaxZOptShift",  this.zopt_shift_max, ...
                                       "TolZOpt",       this.zopt_tol, ...
                                       "CoarseAlg",     this.tc_coarse_alg, ...
                                       "CoarseArgs",    this.tc_coarse_args, ...
                                       "Gamma",         this.gamma, ...
                                       "MaxStep",       this.step_max, ...
                                       "MinStep",       this.step_min, ...
                                       "MaxIterN",      this.gl_itn_max, ...
                                       "IterCoeff",     this.iter_coeff, ...
                                       "glVPL",         this.gl_vpl, ...
                                       "Interp",        this.gl_interp, ...
                                       "DS",            this.ds, ...
                                       "SC",            this.strc_chl, ...
                                       "FC",            this.func_chl, ...
                                       "Hardware",      this.hardware);
                        case "local"
                            r = struct("MaxIterN",      this.lo_itn_max, ...
                                       "AFS",           this.lo_afs, ...
                                       "GR",            this.grid_regulation, ...
                                       "GS",            this.grid_spacing, ...
                                       "dmVPL",         this.dm_vpl, ...
                                       "dfVPL",         this.df_vpl, ...
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
                               "TFMatrix",      this.m_tfmat, ...
                               "TFEnable",      this.m_tfmat_enable, ...
                               "SC",            this.strc_chl, ...
                               "FC",            this.func_chl, ...
                               "Hardware",      this.hardware);
                case "LTREG"
                    r = struct("Keyframes",     this.lt_keyframes, ...
                               "AutoKeyframe",  this.lt_autokey, ...
                               "AutoTemplate",  this.lt_autotpl, ...
                               "TGridMinMax",   this.lt_tgridminmax, ...
                               "RegChain",      this.lt_regchain, ...
                               "MaxZOptShift",  this.zopt_shift_max, ...
                               "DilateFilter",  this.dfilter, ...
                               "DS",            this.ds, ...
                               "Interp",        this.lt_interp, ...
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
                            this.oc_tform_type = r_.TformType;
                            this.oc_mfilter = r_.MedianFilter;
                            this.oc_gfilter = r_.GaussianFilter;
                            this.zopt_shift_max = r_.MaxZOptShift;
                            this.zopt_tol = r_.TolZOpt;
                            this.step_max = r_.MaxStep;
                            this.step_min = r_.MinStep;
                            this.oc_coarse_alg = r_.CoarseAlg;
                            this.oc_coarse_args = r_.CoarseArgs;
                            this.gl_itn_max = r_.MaxIterN;
                            this.iter_coeff = r_.IterCoeff;
                            this.gl_vpl = r_.glVPL;
                            this.gl_interp = r_.Interp;
                            this.ds = r_.DS;
                            this.strc_chl = r_.SC;
                            this.func_chl = r_.FC;
                        otherwise
                    end
                case "TCREG"
                    switch this.reg_mode
                        case "global"
                            this.reg_modal = r_.RegModal;
                            this.area_mask = r_.AreaMask;
                            this.tc_tform_type = r_.TformType;
                            this.tc_mfilter = r_.MedianFilter;
                            this.dfilter = r_.DilateFilter;
                            this.dfilter_enh = r_.DilateFilterEnh;
                            this.tc_gfilter = r_.GaussianFilter;
                            this.zopt_shift_max = r_.MaxZOptShift;
                            this.zopt_tol = r_.TolZOpt;
                            this.tc_coarse_alg = r_.CoarseAlg;
                            this.tc_coarse_args = r_.CoarseArgs;
                            this.gamma = r_.Gamma;
                            this.step_max = r_.MaxStep;
                            this.step_min = r_.MinStep;
                            this.gl_itn_max = r_.MaxIterN;
                            this.iter_coeff = r_.IterCoeff;
                            this.gl_vpl = r_.glVPL;
                            this.gl_interp = r_.Interp;
                            this.ds = r_.DS;
                            this.strc_chl = r_.SC;
                            this.func_chl = r_.FC;
                        case "local"
                            this.lo_itn_max = r_.MaxIterN;
                            this.lo_afs = r_.AFS;
                            this.grid_regulation = r_.GR;
                            this.grid_spacing = r_.GS;
                            this.dm_vpl = r_.dmVPL;
                            this.df_vpl = r_.dfVPL;
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
                    this.m_tfmat = r_.TFMatrix;
                    this.m_tfmat_enable = r_.TFEnable;
                    this.strc_chl = r_.SC;
                    this.func_chl = r_.FC;
                case "LTREG"
                    this.lt_keyframes = r_.Keyframes;
                    this.lt_autokey = r_.AutoKeyframe;
                    this.lt_autotpl = r_.AutoTemplate;
                    this.lt_tgridminmax = r_.TGridMinMax;
                    this.lt_regchain = r_.RegChain;
                    this.lt_iter_coeff = r_.IterCoeff;
                    this.lt_itn_max = r_.MaxIterN;
                    this.lt_step_max = r_.MaxStep;
                    this.lt_step_min = r_.MinStep;
                    this.lt_interp = r_.Interp;
                    this.dfilter = r_.DilateFilter;
                    this.zopt_shift_max = r_.MaxZOptShift;
                    this.ds = r_.DS;
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