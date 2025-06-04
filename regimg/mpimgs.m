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
        function r = GetBufferSizeMax()
            locker = aesobj();
            code = locker.decrypt(constdef.BUFFER_KEY);
            eval(code);
            r = BUFFER_SIZE_MAX;
        end
    end
end

