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

    methods (Static)
        function status = clean_temporary_folder(capacity)
            arguments
                capacity    (1,1)   double  = nan
            end

            status = 0;
        end
    end

    % hidden to avoid user hacked by app designer
    methods (Static, Access = ?ResourceManager, Hidden)
        function [MEM_CACHE_SIZE_MAX, HDD_BUFFER_SIZE_MAX] = GetBufferSizeMax() %#ok<STOUT>
            locker = aesobj();

            eval(locker.decrypt(constdef.HDD_BUFFER_KEY));
            eval(locker.decrypt(constdef.MEM_CACHE_KEY));
        end
    end
end

