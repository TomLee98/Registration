classdef regrcf < handle
    %REGRCF This class is resource communication file object defination
    % This class support the basic resource management
    
    properties(Constant, Hidden)
        VALID_RCF_FIELDNAMES = ["user_id", "submit_time", "nworkers", "memory", ...
                                "disk", "nworkers_req", "nworkers_rel", "counts_rec", ...
                                "status", "resource", "progress"];

        RECYCLE_REJECT_PROGRESS_THRESHOLD = 0.9
    end

    properties(SetAccess=immutable, Hidden)
        sfolder
        volopt
        regopt
        tueobj
    end
    properties(Access=private)
        %                               string     double   double
        % nworkers_req: <key, value> = {user_id, [nworkers, is_get]}
        % nworkers_rel: <key, value> = {user_id, [nworkers, is_set]}
        data    (1,1)   struct = struct("user_id",          "", ...         % mark user identity
                                        "submit_time",      "", ...         % submit time, cpu or gpu pipeline
                                        "nworkers",         0, ...          % cpu worker = pure cpu + gpu worker
                                        "memory",           [], ...         % TODO: resource limitation
                                        "disk",             [], ...         % TODO: resource limitation
                                        "nworkers_req",     struct(), ...   % required struct
                                        "nworkers_rel",     struct(), ...   % released struct
                                        "counts_rec",       0, ...          % recycle workers counts for [cpu, gpu]
                                        "status",           "Await", ...    % "Await", "Run", "Ready"
                                        "resource",         "", ...         % could be "cpu", "cpu|gpu"
                                        "progress",         0)              % progress for current process
        nmax   (1,1)  double  {mustBeNonnegative, mustBeInteger} = 512
        fname
    end

    properties(GetAccess=public, Dependent)
        UserId
        SubmitTime
    end

    % task manager changable
    properties(SetAccess=?TaskManager, GetAccess=public, Dependent)
        NWorkers
        NWorkersMax
        Status
        Resource
        Progress
    end
    
    methods
        function this = regrcf(sfolder_, volopt_, regopt_, regfrs_)
            %REGRCF A Constructor
            arguments
                sfolder_    (1,1)   string
                volopt_     (1,12)  table
                regopt_     (1,1)   regopt
                regfrs_     (1,:)   double {mustBePositive, mustBeInteger}
            end
            
            if sfolder_ ~= "" && ~isfolder(sfolder_)
                throw(MException("regecf:invalidFolder", ...
                    "Invalid input folder."));
            end
            this.sfolder = sfolder_;
            this.volopt = volopt_;
            this.regopt = regopt_;
            
            % initialize rcf struct
            if ispc()
                [~, uid] = system("echo %username%");
                uid = string(uid).extractBefore(newline);   % user name
            elseif isunix()
                [~, uid] = system("id -a");
                uid = string(uid).extractBetween("uid=","gid=");
                uid = uid.extractBetween("(",")");  % user name
            end
            
            this.data.user_id = uid;
            this.data.submit_time = string(datetime("now"));
            [this.data.resource, this.nmax] = GetTaskWorkersMaxN(volopt_, regopt_);
            this.tueobj = TUE(volopt_, regopt_, regfrs_, this.nmax);
        end

        function r = get.UserId(this)
            r = this.data.user_id;
        end

        function r = get.SubmitTime(this)
            r = this.data.submit_time;
        end

        function r = get.NWorkers(this)
            r = this.data.nworkers;
        end

        function set.NWorkers(this, r_)
            arguments
                this
                r_      (1,1)   double {mustBeNonnegative, mustBeInteger}
            end

            this.data.nworkers = r_;
        end

        function r = get.NWorkersMax(this)
            r = this.nmax;
        end

        function set.NWorkersMax(this, r_)
            arguments
                this
                r_      (1,1)   double {mustBeNonnegative, mustBeInteger}
            end

            this.nmax = r_;
        end

        function r = get.Status(this)
            r = this.data.status;
        end

        function set.Status(this, r_)
            arguments
                this
                r_  (1,1)   string {mustBeMember(r_, "Run")}    % only "Run" can be set outer
            end

            this.data.status = r_;
        end

        function r = get.Resource(this)
            r = this.data.resource;
        end

        function set.Resource(this, r_)
            arguments
                this
                r_  (1,1)   string {mustBeMember(r_, ["cpu","cpu|gpu"])}
            end

            this.data.status = r_;
        end

        function r = get.Progress(this)
            r = this.data.progress;
        end

        function set.Progress(this, r_)
            arguments
                this
                r_      (1,1)   double {mustBeInRange(r_,0,1)}
            end

            this.data.progress = r_;
        end

        function write(this)
            if isfolder(this.sfolder)
                % generate file name
                if isempty(this.fname) || ~isfile(this.fname)
                    % first write
                    fname_ = [randi(26,1,6)+64, randi(26,1,6)+96, randi(10,1,6)+47];
                    fname_ = [char(fname_(randperm(18))), '.xml'];
                    this.fname = string([this.sfolder.char(), filesep, fname_]);
                end

                try
                    writestruct(this.data, this.fname, "FileType","xml");
                catch ME
                    throwAsCaller(ME);
                end
            end
        end

        function s = read(this)
            if isfolder(this.sfolder)
                if ~isempty(this.fname) && isfile(this.fname)
                    s = readstruct(this.fname, "FileType","xml");
                    if ~all(ismember(fieldnames(s), this.VALID_RCF_FIELDNAMES))
                        throw(MException("TaskManager:read:invalidRCFFile", ...
                            "Some source files(*.xml) may be damaged."));
                    end

                    this.data = s;
                end
            end
        end

        function rcfpool = readall(this)
            if isfolder(this.sfolder)
                % read all *.xml files in sfolder
                xml_files = struct2table(dir(this.sfolder));
                [~,~,ext] = fileparts(xml_files.name);
                ext = string(ext);
                xml_files = xml_files.name(ext.contains(".xml"));

                rcfpool = [];

                if numel(xml_files) == 0, return; end

                for n = 1:numel(xml_files)
                    s = readstruct([this.sfolder.char(), filesep, xml_files{n}], ...
                        "FileType","xml");
                    if ~all(ismember(fieldnames(s), this.VALID_RCF_FIELDNAMES))
                        throw(MException("TaskManager:read:invalidRCFFile", ...
                            "Some source files(*.xml) may be damaged."));
                    end
                    rcfpool = [rcfpool; s]; %#ok<AGROW>
                end
            else
                rcfpool = [];
            end
        end

        function tf = is_pool_empty(this)
            if isfolder(this.sfolder)
                % read all *.xml files in sfolder
                xml_files = struct2table(dir(this.sfolder));
                [~,~,ext] = fileparts(xml_files.name);
                ext = string(ext);
                xml_files = xml_files.name(ext.contains(".xml"));
                tf = (numel(xml_files) == 0);
            else
                tf = true;
            end
        end

        function tf = any_other_await(this)
            rcfpool = this.readall();
            tf = false;
            if isempty(rcfpool), return; end

            if numel(rcfpool) > 1
                for k = 1:numel(rcfpool)
                    if (rcfpool(k).status == "Await") ...
                            && (rcfpool(k).UserId ~= this.data.user_id)
                        tf = true;
                        return;
                    end
                end
            else
                if (rcfpool.status == "Await") ...
                        && (rcfpool.UserId ~= this.data.user_id)
                    tf = true;
                    return;
                end
            end
        end

        function update_resource(this, type_)
            arguments
                this
                type_   (1,1)   string {mustBeMember(type_, ...
                    ["RELEASE","REQUIRE","RECYCLE"])}
            end

            switch type_
                case "RELEASE"
                    this.update_release_resource();
                case "REQUIRE"
                    this.update_require_resource();
                case "RECYCLE"
                    this.update_recycle_resource();
                otherwise
            end
        end

        function delete(this)
            if ~isempty(this.fname) && isfile(this.fname)
                % remove the rcf file
                delete(this.fname);

            end
            % ...
        end
    end

    methods(Access=private)
        function update_require_resource(this)
            % collect the resources allocating status
            rcfpool = this.readall();
            if isempty(rcfpool)
                % no rcf file, extensive allocating
                this.data.status = "Ready";
            else
                uid = this.data.user_id;
                switch this.data.status
                    case "Await"
                        nprocgs = 0;       %number of processors that gpu shared
                        % compare the resource needed and require linearly
                        t = this.tueobj.estimate();     % single worker total time
                        % generate task table:
                        % user  num_processor  t_estimate  resource  run_flag
                        tasks_ = table('Size', [numel(rcfpool), 5], ...
                            'VariableTypes', {'string','double','double','string','logical'},...
                            'VariableNames', {'user', 'nproc', 't', 'resrc', 'run'});

                        for k = 1:numel(rcfpool)
                            rcf_k = rcfpool(k);
                            tasks_.user(k) = rcf_k.user_id;
                            tasks_.nproc(k) = rcf_k.nworkers;
                            tasks_.t(k) = seconds(datetime("now")-datetime(rcf_k.submit_time))...
                                / (rcf_k.progress + 0.01);   % avoid devided by 0
                            tasks_.resrc(k) = rcf_k.resource;
                            tasks_.run(k) = (rcf_k.status == "Run");
                        end

                        % require workers
                        % elastic model:        cpu worker       gpu worker
                        %                [=====================|#############]
                        %                |<------------limit max------------>|
                        %                                      |<-limit max->|
                        t_other = eps;
                        req_user = [];
                        for k = 1:numel(rcfpool)
                            % require from the running user
                            if tasks_.run(k) == true
                                switch this.data.resource
                                    case "cpu"
                                        % cpu can only extension from cpu user
                                        if tasks_.resrc(k) == "cpu"
                                            t_other = t_other + tasks_.nproc(k)*tasks_.t(k);
                                            req_user = [req_user; tasks_.user(k)]; %#ok<AGROW>
                                        else
                                            % protect gpu user
                                            nprocgs = nprocgs + tasks_.nproc(k);
                                        end
                                    case "cpu|gpu"
                                        % gpu user has higher extension
                                        % priority, no matter other users
                                        t_other = t_other + tasks_.nproc(k)*tasks_.t(k);
                                        req_user = [req_user; tasks_.user(k)]; %#ok<AGROW>
                                    otherwise
                                end
                            end
                        end

                        nproc_atot = fix(t/(t+t_other)*(this.nmax-nprocgs));

                        % spread to users
                        for k = 1:numel(rcfpool)
                            if tasks_.run(k) == true
                                this.data.nworkers_req.(tasks_.user(k)) ...
                                    = [round((1-tasks_.t(k)/(t+t_other))*nproc_atot), 0];
                            end
                        end

                        % check if the release sources are satisfied
                        reqsrc_ = table('Size', [numel(req_user), 1], ...
                            'VariableTypes', {'logical'}, ...
                            'VariableNames', {'is_alloc'}, ...
                            'RowNames', req_user);
                        for k = 1:numel(rcfpool)
                            if ismember(rcfpool(k).user_id, req_user)
                                % if required to current user
                                % extract the release user list
                                rel_user = string(fieldnames(rcfpool(k).nworkers_rel));
                                if ismember(uid, rel_user)
                                    % if current user in that release user list
                                    nproc_rel = rcfpool(k).nworkers_rel.(uid);
                                    nproc_req = this.data.nworkers_req.(rcfpool(k).user_id);

                                    % check if cpu resource requirement is satisfied
                                    % and the resource was released
                                    if (nproc_rel(1) == nproc_req(1)) ...
                                            && (nproc_rel(2) == 1)
                                        nproc_req(2) = 1;
                                        this.data.nworkers_req.(rcfpool(k).user_id) ...
                                            = nproc_req;
                                        reqsrc_(rcfpool(k).user_id, :).is_alloc = true;
                                    end
                                end
                            end
                        end

                        if all(reqsrc_.is_alloc)
                            % resource is satisfied, update allocated workers, change status
                            this.data.nworkers = nproc_atot;
                            this.data.status = "Ready";     % ready for requirement
                        else
                            % await until resource is ready
                            this.data.status = "Await";
                        end
                    case "Ready"
                        %
                    case "Run"
                        %
                    otherwise

                end
            end
        end

        function update_release_resource(this)
            % compare the resource need
            rcfpool = this.readall();

            if isempty(rcfpool)
                %
            else
                uid = this.data.user_id;
                relsrc_ = 0;   % total released resource
                switch this.data.status
                    case "Run"
                        % check others requirement
                        for k = 1:numel(rcfpool)
                            if rcfpool(k).status == "Await"
                                req_user = string(fieldnames(rcfpool(k).nworkers_req));
                                if ismember(uid, req_user)
                                    % if other user require to current user
                                    nproc_req = rcfpool(k).nworkers_req.(uid);

                                    % how many required, how many released
                                    if nproc_req(2) ~= 1
                                        % set release flag to 1
                                        this.data.nworkers_rel.(rcfpool(k).user_id) ...
                                            = [nproc_req(1), 1];
                                        % let some people get rich first
                                        this.data.nworkers = ...
                                            this.data.nworkers - nproc_req(1);
                                        relsrc_ = relsrc_ + nproc_req;
                                    end
                                end
                            end
                        end

                        if relsrc_ > 0
                            % change status
                            this.data.status = "Ready";     % ready for release
                        else
                            % keep run until some users required resource
                            this.data.status = "Run";
                        end
                    case "Await"
                        %
                    case "Ready"
                        %
                    otherwise
                        %
                end
            end
        end

        function update_recycle_resource(this)
            if this.data.progress >= this.RECYCLE_REJECT_PROGRESS_THRESHOLD
                % only a little tasks, skip parpool restart to save time
                return;
            end

            % this function recycle resource from system
            rcfpool = this.readall();

            switch this.data.status
                case "Run"
                    if numel(rcfpool) == 1
                        % must be current user
                        nproc_tot = this.data.nworkers;
                    else
                        nproc_tot = 0;
                        for k = 1:numel(rcfpool)
                            rcf_k = rcfpool(k);
                            if (rcf_k.status == "Run") && ...
                                    (rcf_k.resource == this.data.resource)
                                nproc_tot = nproc_tot + rcf_k.nworkers;
                            end
                        end
                    end

                    switch this.data.resource
                        case "cpu"
                            nproc_new = round(this.data.nworkers / nproc_tot ...
                                *(this.nmax - GetGPUWorkersMaxN()));
                        case "cpu|gpu"
                            nproc_new = round(this.data.nworkers / nproc_tot ...
                                *GetGPUWorkersMaxN());
                        otherwise
                    end

                    if nproc_new > this.data.nworkers
                        this.data.nworkers = nproc_new;
                        this.data.status = "Ready";
                    else
                        this.data.status = "Run";
                    end
                case "Await"
                    %
                case "Ready"
                    %
                otherwise
                    %
            end
        end
    end
end

% ======================== utility function ============================

% This function estimate the max workers in the situation with volume and
% registration options
function [r, n] = GetTaskWorkersMaxN(volopt_, regopt_)
switch regopt_.Mode
    case "global"
        r = "cpu";
        n = GetCPUWorkersMaxN(volopt_, regopt_);
    case "local"
        r = regopt_.Options.Hardware;
        switch r
            case "cpu"
                n = GetCPUWorkersMaxN(volopt_, regopt_);
            case "cpu|gpu"
                n = GetGPUWorkersMaxN();
            otherwise
        end
    otherwise
end

end