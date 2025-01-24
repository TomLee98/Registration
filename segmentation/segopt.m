classdef segopt
    %SEGOPT This class defined as the segmentor options type
    
    properties
        Property1
    end
    
    methods
        function obj = segopt(segalg_, runmode_)
            %SEGOPT A constructor
            arguments
                segalg_     (1,1)   string  {mustBeMember(segalg_, ["MP", "ML"])} = "MP"
                runmode_    (1,1)   string  {mustBeMember(runmode_, ["train", "predict"])} = "predict"
            end

            
        end
    end
end

