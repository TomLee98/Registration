classdef taskParser < handle
    %TASKPARSER This class defines taskParser object, which could map the 
    % (volopt, regopt, regfrs) -> task queue
    % As example, long term rregistration need three steps, and taskParser
    % will generate the pipeline automatically

    properties(Constant, Hidden)
        % shared with TUE
        C_DISTRIBUTION = 2;
        C_EXCLUSIVE = 4;
    end
    
    properties(GetAccess=public, SetAccess=private)
        Results     (1,:)   mQueue = mQueue()
    end

    properties(Access=private, Hidden)
        movtmpl
        regfrs
        volopt
        regopt
        nworker
        distrib
        pfunc
    end
    
    methods
        function this = taskParser(movtmpl_, regfrs_, volopt_, regopt_, blocksz_, distrib_)
            %TASKPARSER A Constructor
            arguments
                movtmpl_    (1,1)   regtmpl 
                regfrs_     (1,:)   double  {mustBePositive, mustBeInteger}
                volopt_     (1,12)  table
                regopt_     (1,1)   regopt
                blocksz_    (1,1)   double  {mustBePositive, mustBeInteger}
                distrib_    (1,1)   logical = false
            end

            this.movtmpl = movtmpl_;
            this.regfrs = regfrs_;
            this.volopt = volopt_;
            this.regopt = regopt_;
            this.nworker = blocksz_;
            this.distrib = distrib_;

            %PARSE This function parse the task and generate task queue
            pathToParser = fileparts(mfilename("fullpath"));
            rf = importdata([pathToParser, filesep, 'configuration.ini']);
            rf = string(rf).split(":");
            if ~ismember(this.regopt.Algorithm, rf(:,1))
                throw(MException("taskParser:parse:invalidRegistrationAlgorithm", ...
                        "Can not parse not registered algorithm."));
            end

            this.pfunc = str2func(rf(this.regopt.Algorithm==rf(:,1), 2));
        end
        
        function parse(this)
            if this.distrib == false
                % in exclusive mode, user can adjust the busy index
                batch_sz = this.nworker*this.C_EXCLUSIVE;
            else
                % distribution mode must keep busy index as 2~3 for better
                % runnning efficiency and response sensitivity
                batch_sz = this.nworker*this.C_DISTRIBUTION;
            end

            this.Results = this.pfunc(this.movtmpl, ...
                                      this.regfrs, ...
                                      this.volopt, ...
                                      this.regopt, ...
                                      batch_sz);
        end
    end
end