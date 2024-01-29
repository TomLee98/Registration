classdef Task
    %TASK This class is a value class, which contains registration task
    % information
    
    properties(SetAccess=private, GetAccess=public)
        RegFrames
        RegTemplate
        RegAlgorithm
        RegOptions
        Status
    end
    
    methods
        function this = Task(regfrs_, tmplfr_, regalg_, regopt_, status_)
            %TASK A Constructor
            arguments
                regfrs_ (1,:) double {mustBePositive, mustBeInteger}
                tmplfr_ (1,1) struct
                regalg_ (1,1) string {mustBeMember(regalg_, ...
                    ["TCREG", "OCREG", "MANREG", "LTREG"])}
                regopt_ (1,1) struct
                status_ (1,1) string {mustBeMember(status_, ...
                    ["Await","Run","Done"])} = "Await"
            end

            % double check arguments
            if (numel(unique(regfrs_)) ~= numel(regfrs_)) ...
                    || any(diff(regfrs_) <= 0)
                warning("Task:invalidRegistrationFrames", ...
                    "Task was refined.");
                regfrs_ = sort(unique(regfrs_));
            end

            if ~isempty(setxor(fieldnames(tmplfr_), {'L','G'})) ...
                    || (numel(tmplfr_.L)~=2) || (numel(tmplfr_.G)~=2)
                throw(MException("Task:invalidTemplate", ...
                    "Template must contains global and local information."));
            end

            this.RegFrames = regfrs_;
            this.RegTemplate = tmplfr_;
            this.RegAlgorithm = regalg_;
            this.RegOptions = regopt_;
            this.Status = status_;
        end
    end
end

