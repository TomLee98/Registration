classdef Estimator < handle
    %SIGNALEXTRACTOR This class define a signal extractor object

    properties(Constant, Hidden)
    end

    properties(GetAccess=public, Dependent)
        Activities

    end
    
    properties(Access=private, Hidden)
        options         % 1-by-1 sigopt object
        channel         % 1-by-1 string, "r", "g" or "b"
        components      % 1-by-n nonnegtive integer array, component indices
        activities      
    end

    properties(SetAccess=immutable, Hidden)
        image_src      % 1-by-1 regmov object, the calcium image data 
    end
    
    methods
        function this = Estimator(src_, opts_)
            arguments
                src_     (1,1)   regmov
                opts_    (1,1)   sigopt
            end

            this.image_src = src_;
            this.options = opts_;
        end
        
        function fit(this, opts_, ch_, comps_)
            % This function fit the signal model
            arguments
                this
                opts_   (1,1)   sigopt
                ch_     (1,1)   string  {mustBeMember(ch_, ["r","g","b"])}
                comps_  (1,:)   double  {mustBeNonnegative, mustBeInteger}
            end

            this.options = opts_;
            this.channel = ch_;
            this.components = comps_;

            % 
        end

    end
end

