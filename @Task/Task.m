classdef Task
    %TASK This class is a value class, which contains registration task
    % information that can be read by RegisterWorker object
    
    properties(SetAccess=private, GetAccess=public)
        RegFrames
        RegTemplate
        RegOptions
    end

    properties(Access=private)
        status  (1,1) string {mustBeMember(status, ["Await","Run","Done"])} = "Await"
    end

    properties(SetAccess={?TaskManager, ?RegisterWorker}, GetAccess=?TaskManager, Dependent)
        Status
    end
    
    methods
        function this = Task(regfrs_, tmplfr_, regopt_)
            %TASK A Constructor
            arguments
                regfrs_ (1,:) double {mustBePositive, mustBeInteger}
                tmplfr_ (1,1) regtmpl
                regopt_ (1,1) regopt
            end

            % double check arguments
            if (numel(unique(regfrs_)) ~= numel(regfrs_)) ...
                    || any(diff(regfrs_) <= 0)
                warning("Task:invalidRegistrationFrames", ...
                    "Task was refined.");
                regfrs_ = sort(unique(regfrs_));
            end

            this.RegFrames = regfrs_;
            this.RegTemplate = tmplfr_;
            this.RegOptions = regopt_;
        end

        function this = set.Status(this, r_)
            arguments
                this
                r_ (1,1) string {mustBeMember(r_, ...
                    ["Await","Run","Done"])} = "Await"
            end
            
            this.status = r_;
        end

        function r = get.Status(this)
            r = this.status;
        end
    end
end

