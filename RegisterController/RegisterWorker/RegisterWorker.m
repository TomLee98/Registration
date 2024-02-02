classdef RegisterWorker < handle
    %REGWORKER This class defines rigister worker, which accept the a Task
    % comes from TaskManager object and finish it
    
    properties(Access = private, Hidden)
        mov_raw
        mov_aligned
    end
    
    methods
        function this = RegisterWorker(mov_raw_, mov_aligned_)
            %REGWORKER A Constructor
            arguments
                mov_raw_        (1,1)   regmov
                mov_aligned_    (1,1)   regmov
            end

            this.mov_raw = mov_raw_;
            this.mov_aligned = mov_aligned_;
        end
        
        function status = correct(this, task_)
            %This function parse task_ and do registration
            
        end
    end
end

