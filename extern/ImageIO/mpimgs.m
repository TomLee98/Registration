classdef mpimgs < handle
    %MPIMGS This class is mpimg container, which offers a big image file
    % management solution in low memory mechine
    
    properties(Access=private, Hidden)
        mpimg_workers      
    end
    
    methods
        function this = mpimgs()
            %MPIMGS A constructor

        end
        
    end
end

