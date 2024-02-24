classdef TUE < handle
    %TUE This class is Time Used Estimater

    properties(Constant, Hidden)
        % Global REGISTRATION HIDDEN PARAMETERS
        IMREGTFORM_TRS_PPS = 2.60E5                % pixels per second
        IMWARP_TRS_CUBIC_PPS = 6.52E5
        IMWARP_TRS_LINEAR_CUBIC_RATIO = [3, 1];
        IMREGTFORM_TRS_RID_AFF_RATIO = [3,2,1]     % processing speed ratio: imregtform

        % LOCAL REGISTRATION HIDDEN PARAMETERS
        IMWARP_DEMON_TUSE = 10                      % typical using, actually can be omitted
        IMREGDEMONS_ACF_1_PPS = 1.04E5              % pixels per second

        % IMAGE_SIZE CONSTANT, SAME WITH <RESAMPLE>
        EC50 = 256
    end

    properties(SetAccess=immutable, Hidden)
        volopt
        regopt
        nregfr
    end
    
    methods
        function this = TUE(volopt_, regopt_, regfrs_)
            %TUE A Constructor
            arguments
                volopt_     (1,12)  table
                regopt_     (1,1)   regopt
                regfrs_     (1,:)   double {mustBePositive, mustBeInteger}
            end

            this.volopt = volopt_;
            this.regopt = regopt_;
            this.nregfr = numel(regfrs_);
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
                        if isMATLABReleaseOlderThan("R2022b")
                            t = this.imregopzr_time_use() ...
                                + this.imregtform_time_use() ...
                                + 2*this.imwarp_time_use();
                        else
                            t = this.imregmoment_time_use() ...
                                + this.imregtform_time_use() ...
                                + 2*this.imwarp_time_use();
                        end
                    elseif this.regopt.Mode == "local"
                        if this.volopt.Hardware == "cpu"

                        else
                            
                        end
                    end
                case "MANREG"

                case "LTREG"
                    if this.volopt.Hardware == "cpu"

                    else

                    end
                otherwise
            end
        end
    end

    methods(Access=private)

        function t_use = imregtform_time_use(this)
            switch this.regopt.RegType
                case "auto"
                    if this.regopt.Mode == "Global"
                        [~, rp] = ismember(this.regopt.TformType, ...
                            ["translation","rigid","affine"]);
                        procsp = this.IMREGTFORM_TRS_PPS ...
                            *this.IMREGTFORM_TRS_RID_AFF_RATIO(rp);
                        t_use = this.pixels_calc_num() / procsp;
                    else
                        t_use = 0;
                    end
                case "longterm"

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

        function t_use = imregmoment_time_use(this)

        end

        function t_use = imregdeform_time_use(this)
            
        end

        function t_use = imregopzr_time_use(this)

        end

        function t_use = imwarp_time_use(this)
            [~, rp_type] = ismember(this.regopt.TformType, ...
                ["translation","rigid","affine"]);

            switch this.regopt.Mode
                case "Global"
                    [~, rp_alg] = ismember(this.regopt.Interp_Rigid, ["linear","cubic"]);
                    procsp = this.IMWARP_TRS_CUBIC_PPS ...
                        *this.IMWARP_TRS_LINEAR_CUBIC_RATIO(rp_alg) ...
                        *this.IMREGTFORM_TRS_RID_AFF_RATIO(rp_type);
                    t_use = this.pixels_calc_num() / procsp;
                case "Local"
                    t_use = this.IMWARP_DEMON_TUSE;
            end
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

        end
    end
end

