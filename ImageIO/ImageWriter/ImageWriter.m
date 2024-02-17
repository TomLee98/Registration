classdef ImageWriter < handle
    %IMAGEWRITER This class is an image writer defination, which  support
    %mutiple image format writting
    
    properties(Access=private, Hidden)
        caller
        file
    end
    
    methods
        function this = ImageWriter(caller_, file_)
            arguments
                caller_ (1,1)   Register
                file_   (1,1)   string
            end

            this.caller = caller_;
            this.file = file_;
        end
        
        function outputArg = save(this, mov, metadata, block)
            arguments
                this
                mov         (1,1)   regmov
                metadata    (1,1)   struct
                block       (1,1)   double {mustBePositive, mustBeInteger} = 50
            end
            
            
        end
    end
end

