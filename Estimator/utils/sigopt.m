classdef sigopt
    %SIGOPT This class is signal extractor option defination, which is packaged
    %value class
    
    properties(Access=private, Hidden)
        bl_model    (1,1) string {mustBeMember(bl_model, ["MovingQuantile", "MixedExponential"])} = "MovingQuantile"
        bkg_model   (1,1) string {mustBeMember(bkg_model, ["Constant", "Flexible"])} = "Constant"
        q_mm_min    (1,1) double {mustBeInRange(q_mm_min, 0, 100)} = 10
        q_me_min    (1,1) double {mustBeInRange(q_me_min, 0, 100)} = 10
        q_mm_auto   (1,1) logical = false
        win_size    (1,1) double {mustBePositive, mustBeInteger} = 100
        order       (1,1) double {mustBeInteger, mustBeInRange(order, 1,5)} = 1
        order_auto  (1,1) logical = false
        background  (1,1) double {mustBeNonnegative} = 100
        kernel      (1,1) string {mustBeMember(kernel, ["uniform","gaussian","log"])} = "uniform"
        kernel_auto (1,1) logical = false
        f_only      (1,1) logical = false
        noise_model (1,1) string {mustBeMember(noise_model, ["normal","exponential","gamma","none"])} = "normal"
        sig_level   (1,1) string {mustBeMember(sig_level, ["0.05", "0.01", "0.001", "0.0001"])} = "0.01"
    end

    properties(Access=public, Dependent)
        BaselineModel
        BackgroundModel
        Options
    end
    
    methods
        function this = sigopt(baseline_, background_)
            %SIGOPT A Constructor
            arguments
                baseline_   (1,1)   string {mustBeMember(baseline_, ["MovingQuantile", "MixedExponential"])} = "MovingQuantile"
                background_ (1,1)   string {mustBeMember(background_, ["Constant", "Flexible"])} = "Constant"
            end
            this.bl_model = baseline_;
            this.bkg_model = background_;
        end

        function r = get.BackgroundModel(this)
            r = this.bkg_model;
        end

        function this = set.BackgroundModel(this, r)
            arguments
                this
                r   (1,1)   string {mustBeMember(r, ["Constant", "Flexible"])} = "Constant"
            end

            this.bkg_model = r;
        end

        function r = get.BaselineModel(this)
            r = this.bl_model;
        end

        function this = set.BaselineModel(this, r)
            arguments
                this
                r   (1,1)  string {mustBeMember(r, ["MovingQuantile", "MixedExponential"])} = "MovingQuantile"
            end

            this.bl_model = r;
        end

        function r = get.Options(this)
            r = this.gen_option_();
        end
    end

    methods(Access=public, Hidden)
        function this = set(this, varargin)
            p = inputParser;
            p.StructExpand = false;     % allow parameters struct

            % =========== not overloading parameters =============
            addParameter(p, 'WindowSize',       this.win_size);
            addParameter(p, 'RegressOrder',     this.order);
            addParameter(p, 'AutoRegress',      this.order_auto);
            addParameter(p, 'Background',       this.background);
            addParameter(p, 'Kernel',           this.kernel);
            addParameter(p, 'AutoKernel',       this.kernel_auto);
            addParameter(p, 'FluorescenceOnly', this.f_only);
            addParameter(p, 'NoiseModel',       this.noise_model);
            addParameter(p, 'Significance',     this.sig_level);

            switch this.bl_model
                case "MovingQuantile"
                    addParameter(p, 'MinQuantileValue', this.q_mm_min);
                    addParameter(p, 'AutoQuantile',     this.q_mm_auto);
                case "MixedExponential"
                    addParameter(p, 'MinQuantileValue', this.q_me_min);
                otherwise
            end

            switch this.bkg_model
                case "Constant"
                    
                case "Flexible"

                otherwise
            end

            parse(p, varargin{:});

            this = this.set_option_(p.Results);
        end
    end

     methods(Access=private)
        function r = gen_option_(this)
            r = struct("WindowSize",        this.win_size, ...
                       "RegressOrder",      this.order, ...
                       "AutoRegress",       this.order_auto, ...
                       "Background",        this.background, ...
                       "Kernel",            this.kernel, ...
                       "AutoKernel",        this.kernel_auto, ...
                       "FluorescenceOnly",  this.f_only, ...
                       "NoiseModel",        this.noise_model, ...
                       "Significance",      this.sig_level);

            % dynamic items added
            switch this.bl_model
                case "MovingQuantile"
                    r.MinQuantileValue = this.q_mm_min;
                    r.AutoQuantile = this.q_mm_auto;
                case "MixedExponential"
                    r.MinQuantileValue = this.q_me_min;
                otherwise
            end
        end

        function this = set_option_(this, r_)
            % r_ is parser results struct, inverse the mapping

            switch this.bl_model
                case "MovingQuantile"
                    this.q_mm_min = r_.MinQuantileValue;
                    this.q_mm_auto = r_.AutoQuantile;
                case "MixedExponential"
                    this.q_me_min = r_.MinQuantileValue;
                otherwise
            end

            this.win_size = r_.WindowSize;
            this.order = r_.RegressOrder;
            this.order_auto = r_.AutoRegress;
            this.background = r_.Background;
            this.kernel = r_.Kernel;
            this.kernel_auto = r_.AutoKernel;
            this.f_only = r_.FluorescenceOnly;
            this.noise_model = r_.NoiseModel;
            this.sig_level = r_.Significance;
        end
    end
end

