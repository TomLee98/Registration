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
            arguments
                this
                task_   (1,1)   task
            end

            regfrs_ = task_.RegFrames;
            Algorithm = task_.RegOptions.Algorithm;
            regopt_ =  task_.RegOptions.Options;
            regopt_.Mode = task_.RegOptions.Mode;   % dynamic combine registration mode
            regopt_.SubAlgorithm = task_.RegOptions.SubAlgorithm;
            movtmpl_ = task_.RegTemplate.RefVol;

            switch Algorithm
                case "OCREG"
                    status = ocreg(this.mov_raw, this.mov_aligned, movtmpl_, ...
                                    regfrs_, regopt_);
                case "TCREG"
                    status = tcreg(this.mov_raw, this.mov_aligned, movtmpl_, ...
                                    regfrs_, regopt_);
                case "LTREG"
                    status = ocreg(this.mov_raw, this.mov_aligned, movtmpl_, ...
                                    regfrs_, regopt_);
                case "MANREG"
                    status = manreg(this.mov_raw, this.mov_aligned, movtmpl_, ...
                                    regfrs_, regopt_);
                otherwise
                    throw(MException("RegisterWorker:unregisteredFunction", ...
                        "Unsupported registration algorithm."));
            end
        end

        function delete(this)
            %
        end
    end
end

