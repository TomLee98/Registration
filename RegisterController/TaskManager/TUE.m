classdef TUE < handle
    %TUE This class is Time Used Estimater

    properties(Constant, Hidden)
        % Global REGISTRATION HIDDEN PARAMETERS
        IMREGTFORM_TRANSL_PPS = 1E5       % pixels per second
        IMREGTFORM_RIGID_PPS = 3E4        % pixels per second
        IMREGTFORM_AFFINE_PPS = 2.5E4       % pixels per second

        IMREGDEMONS_GPU_PPS = 1.5E4           % @ AFS = 1.5
        IMREGDEMONS_CPU_PPS = 1.5E3           % TODO: test the real acceleration ratio, may be not 10:1
        IMHISTMATCHN_PPS = 1E5              % @ LOW LEVEL

        IMREGOPZR_PPS = 3.5E5                 % pixels per second

        PREPROC_PPS = 3.5E5                   % pixels per second

        IMWARP_LINEAR_SPO = 0.05            % seconds per operation, nearly O(1) vs pixels number
        IMWARP_CUBIC_SPO = 0.3              % seconds per operation, ...

        % IMAGE_SIZE CONSTANT, SAME WITH <RESAMPLE>
        EC50 = 256

        YIELD_TIME_CONST = 0.1              % parallel toolbox configuration
        PARALLEL_CONSTANT = 7/8             % parallel toolbox configuration

        C_DISTRIBUTION = 2;
        C_EXCLUSIVE = 4;
    end

    properties(SetAccess=immutable, Hidden)
        volopt
        regopt
        nregfr
        nws
    end
    
    methods
        function this = TUE(volopt_, regopt_, regfrs_, nws_)
            %TUE A Constructor
            arguments
                volopt_     (1,12)  table
                regopt_     (1,1)   regopt
                regfrs_     (1,:)   double {mustBePositive, mustBeInteger}
                nws_        (1,1)   double {mustBePositive, mustBeInteger}
            end

            this.volopt = volopt_;
            this.regopt = regopt_;
            this.nregfr = numel(regfrs_);
            this.nws = nws_;
        end
        
        function t = estimate(this, distrib)
            % this function generate the resource level mapping:
            % (volopt, regopt) -> positive real number (normalized calculation time)
            t = 0;
            switch this.regopt.Algorithm
                case "OCREG"
                    if this.regopt.Mode == "global"

                    elseif this.regopt.Mode == "local"
                        if this.regopt.Options.Hardware == "cpu"
                            
                        elseif this.regopt.Options.Hardware == "cpu|gpu"
                            
                        end
                    end
                case "TCREG"
                    if this.regopt.Mode == "global"
                        % in global algorithm, tcreg has pipeline with:
                        % about 1X proproc, 1X imregopzr, 1X imregtform, 
                        % 2X imwarp, and 1X independent yield time estimation, 
                        % so time calculation as:
                        t = this.preproc_time_use() + ...
                            this.imregopzr_time_use() + ...
                            this.imregtform_time_use() + ...
                            2*this.imwarp_time_use() + ...
                            this.yield_time_use(distrib);
                    elseif this.regopt.Mode == "local"
                        % in local algorithm, tcreg has pipeline with:
                        % 1X imregdeform or 1X imregdemons, 1X imhistmatchn(optional), 
                        % 2X imwarp(can be omitted), and 1X independent 
                        % yield time estimation, 
                        % so time calculation as:
                        if this.regopt.SubAlgorithm == "usual"
                            % imregdeform case
                             t = this.imregdeform_time_use() + ...
                                 this.imhistmatchn_time_use();
                        elseif this.regopt.SubAlgorithm == "advanced"
                            % imregdemons case
                             t = this.imregdemons_time_use() + ...
                                 this.imhistmatchn_time_use();
                        end
                    end
                case "MANREG"

                case "LTREG"

                otherwise
            end
        end
    end

    methods(Access=private)

        function t_use = imregtform_time_use(this)
            switch this.regopt.Options.TformType
                case "translation"
                    t_use = this.pnpc()/this.IMREGTFORM_TRANSL_PPS;
                case "rigid"
                    t_use = this.pnpc()/this.IMREGTFORM_RIGID_PPS;
                case "affine"
                    t_use = this.pnpc()/this.IMREGTFORM_AFFINE_PPS;
                otherwise
                    t_use = 0;
            end
        end

        function t_use = imregdemons_time_use(this)
            % TODO: use more accuracy time estimation
            if this.volopt.Hardware == "cpu"
                t_use = this.pnpc()/this.IMREGDEMONS_CPU_PPS;   % @AFS = 1.5
            elseif this.volopt.Hardware == "cpu|gpu"
                t_use = this.pnpc()/this.IMREGDEMONS_GPU_PPS;   % @AFS = 1.5
            end
        end

        function t_use = imregdeform_time_use(this)
            t_use = 0;
        end

        function t_use = imhistmatchn_time_use(this)
            if this.regopt.Options.ImageRehist == true
                t_use = this.pnpc()/this.IMHISTMATCHN_PPS;  % @LOW-LEVEL
            else
                t_use = 0;
            end
        end

        function t_use = imregopzr_time_use(this)
            t_use = this.pnpc()/this.IMREGOPZR_PPS;
        end

        function t_use = imwarp_time_use(this)
            switch this.regopt.Options.Interp
                case "linear"
                    t_use = this.nregfr*this.IMWARP_LINEAR_SPO;
                case "cubic"
                    t_use = this.nregfr*this.IMWARP_CUBIC_SPO;
                otherwise
            end
        end

        function t_use = preproc_time_use(this)
            t_use = this.pnpc()/this.PREPROC_PPS;
        end

        function t_use = yield_time_use(this, distrib)
            % an experience formation:
            % T_N = (N_CPU/C*2)/w*(N+100)*T
            % note that: This is a typical estimation. The real situation is
            % dependents on some frames convergence speed, which could be 
            % prior estimated hardly

            N_CPU = round(feature('numCores')*this.PARALLEL_CONSTANT);

            if distrib == true
                % tasker manager fixed the batch size
                C = this.C_DISTRIBUTION;
            else
                C = this.C_EXCLUSIVE;  % could be modified by user in some versions
            end

            t_use = this.YIELD_TIME_CONST * ...     % seconds per yield loop
                (N_CPU/C*2)/this.nws*(this.nregfr + 100);
        end

        function pn = pnpc(this)
            % This function calculate the pixels number per channel
            if this.volopt.channels == 1

            elseif this.volopt.channels == 2
                if this.regopt.Mode == "global"
                    if this.regopt.Options.DS == "auto"
                        sz_max = max(this.volopt.width, this.volopt.height);
                        scale = 1/(this.EC50/(this.EC50 + sz_max));
                    else
                        scale = str2double(string(this.regopt.Options.DS).extractBefore(2));
                    end
                elseif this.regopt.Mode == "local"
                    scale = 1;
                end

                % single color channel
                pn =  this.volopt.width*this.volopt.height/scale.^2 ... % XY
                    *this.volopt.slices ...                             % Z
                    *this.nregfr;                                       % T
            end

            pn = ceil(pn);
        end
    end
end

