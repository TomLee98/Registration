classdef TUE < handle
    %TUE This class is Time Used Estimater

    properties(Constant, Hidden)
        % Global REGISTRATION HIDDEN PARAMETERS
        IMREGTFORM_TRANSL_PPS = 3.4E5       % pixels per second
        IMREGTFORM_RIGID_PPS = 8.5E4        % pixels per second
        IMREGTFORM_AFFINE_PPS = 7.5E4       % pixels per second

        IMREGOPZR_PPS = 1E6                 % pixels per second

        PREPROC_PPS = 1E6                   % pixels per second

        IMWARP_LINEAR_SPO = 0.05            % seconds per operation
        IMWARP_CUBIC_SPO = 0.3              % seconds per operation

        % IMAGE_SIZE CONSTANT, SAME WITH <RESAMPLE>
        EC50 = 256

        YIELD_TIME_CONST = 0.1              % from parallel toolbox configuration
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
        
        function t = estimate(this)
            % this function generate the resource level mapping:
            % (volopt, regopt) -> positive real number (normalized calculation time)
            t = 0;
            switch this.regopt.Algorithm
                case "OCREG"
                    if this.regopt.Mode == "global"

                    elseif this.regopt.Mode == "local"
                        if this.volopt.Hardware == "cpu"
                            
                        elseif this.volopt.Hardware == "cpu|gpu"
                            
                        end
                    end
                case "TCREG"
                    if this.regopt.Mode == "global"
                        % in global algorithm, tcreg has pipeline with:
                        % about 1X proproc, 1X imregopzr, 1X imregtform, 
                        % 2X imwarp, and 1X indipenmdent yield time estimation, 
                        % so time calculation as:
                        t = this.preproc_time_use() + ...
                            this.imregopzr_time_use() + ...
                            this.imregtform_time_use() + ...
                            2*this.imwarp_time_use() + ...
                            this.yield_time_use();
                    elseif this.regopt.Mode == "local"

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
                    t_use = this.pixels_calc_num()/this.IMREGTFORM_TRANSL_PPS;
                case "rigid"
                    t_use = this.pixels_calc_num()/this.IMREGTFORM_RIGID_PPS;
                case "affine"
                    t_use = this.pixels_calc_num()/this.IMREGTFORM_AFFINE_PPS;
                otherwise
                    t_use = 0;
            end
        end

        function t_use = imregdemons_time_use(this)
            switch this.regopt.Mode
                case "Local"
                    pn =  this.volopt.width*this.volopt.height*this.volopt.slices*this.nregfr;
                    t_use = pn / this.IMREGDEMONS_ACF_1_PPS;
                    t_use = 0.5*t_use *(this.regopt.AFS+1);
                otherwise
                    t_use = 0;
            end
        end

        function t_use = imregdeform_time_use(this)
            
        end

        function t_use = imregopzr_time_use(this)
            t_use = this.pixels_calc_num()/this.IMREGOPZR_PPS;
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
            t_use = this.pixels_calc_num()/this.PREPROC_PPS;
        end

        function t_use = yield_time_use(this)
            % an experience formation
            t_use = this.YIELD_TIME_CONST * ...
                4/this.nws*(this.nregfr + 100);
        end

        function pn = pixels_calc_num(this)
            if this.volopt.channels == 1

            elseif this.volopt.channels == 2
                if this.regopt.DS == "auto"
                    sz_max = max(this.volopt.width, this.volopt.height);
                    scale = 1/(this.EC50/(this.EC50 + sz_max));
                else
                    scale = str2double(string(this.regopt.DS).extractBefore(2));
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

