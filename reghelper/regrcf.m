classdef regrcf < handle
    %REGRCF This class is resource communication file object defination
    % This class support the basic resource management
    
    properties(Constant, Hidden)
        VALID_RCF_FIELDNAMES = ["user_id", "submit_time", "nworkers", "nbatches", ...
                                "memory", "disk", "nworkers_req", "nworkers_rel", ...
                                "counts_rec", "status", "resource", "progress"];

        RECYCLE_REJECT_PROGRESS_THRESHOLD = 0.9
        SPEED_UP_THRESHOLD = 0.25
        USERS_NUMBER_MAX = 4
    end

    properties(SetAccess=immutable, Hidden)
        sfolder
        volopt
        regopt
        distrib
    end

    properties(Access=private)
        %                               string     double   double
        % nworkers_req: <key, value> = {user_id, [nworkers, is_get]}
        % nworkers_rel: <key, value> = {user_id, [nworkers, is_set]}
        data    (1,1)   struct = struct("user_id",          "", ...         % mark user identity
                                        "submit_time",      "", ...         % submit time, cpu or gpu pipeline
                                        "nworkers",         0, ...          % cpu worker = pure cpu + gpu worker
                                        "nbatches",         0, ...          % total batches in the job
                                        "memory",           [], ...         % TODO: resource limitation
                                        "disk",             [], ...         % TODO: resource limitation
                                        "nworkers_req",     struct("null", 0), ...   % required struct
                                        "nworkers_rel",     struct("null", 0), ...   % released struct
                                        "counts_rec",       0, ...          % recycle workers counts for [cpu, gpu]
                                        "status",           "Await", ...    % "Await", "Run", "Ready"
                                        "resource",         "", ...         % could be "cpu", "cpu|gpu"
                                        "progress",         0)              % progress for current process
        nmax    (1,1)  double  {mustBeNonnegative, mustBeInteger} = 512
        fname
        await_counter   (1,1)   double  = 0
    end

    properties(GetAccess=public, Dependent)
        UserId
        SubmitTime
    end

    % task manager changable
    properties(SetAccess=?TaskManager, GetAccess=public, Dependent)
        NWorkers
        NWorkersMax
        NBatches
        Status
        Resource
        Progress
    end
    
    methods
        function this = regrcf(sfolder_, volopt_, regopt_, distrib_)
            %REGRCF A Constructor
            arguments
                sfolder_    (1,1)   string
                volopt_     (1,12)  table
                regopt_     (1,1)   regopt
                distrib_    (1,1)   logical
            end
            
            if sfolder_ ~= "" && ~isfolder(sfolder_)
                throw(MException("regecf:invalidFolder", ...
                    "Invalid input folder."));
            end
            this.sfolder = sfolder_;
            this.volopt = volopt_;
            this.regopt = regopt_;
            this.distrib = distrib_;
            
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

        function r = get.NBatches(this)
            r = this.data.nbatches;
        end

        function set.NBatches(this, r_)
            arguments
                this
                r_      (1,1)   double {mustBePositive, mustBeInteger}
            end

            this.data.nbatches = r_;
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
                    rng('shuffle');     % make sure a random generation
                    fname_ = [randi(26,1,6)+64, randi(26,1,6)+96, randi(10,1,6)+47];
                    fname_ = [char(fname_(randperm(18))), '.xml'];
                    this.fname = string([this.sfolder.char(), filesep, fname_]);
                end

                try
                    % note that file permission binding on pool folder (groups)
                    % if user without permission, add to group 'regusers'
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
            % This function read all rcf files in pool
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
            % This function determines if current rcf pool is empty
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
            % This function determines if any other user is await, no
            % matter current user status
            rcfpool = this.readall();
            tf = false;
            if isempty(rcfpool), return; end

            if numel(rcfpool) > 1
                for k = 1:numel(rcfpool)
                    if (rcfpool(k).status == "Await") ...
                            && (rcfpool(k).user_id ~= this.data.user_id)
                        tf = true;
                        return;
                    end
                end
            else
                if (rcfpool.status == "Await") ...
                        && (rcfpool.user_id ~= this.data.user_id)
                    tf = true;
                    return;
                end
            end
        end

        function tf = current_await_only(this)
            % This function determines if and only if current user await only,
            % which means the first user join the multi-user scheduling network
            rcfpool = this.readall();
            tf = false;
            if isempty(rcfpool), return; end

            if isscalar(rcfpool)
                if (rcfpool.status == "Await") ...
                        && (rcfpool.user_id == this.data.user_id)
                    tf = true;
                    return;
                end
            end
        end

        function tf = current_await_with_rcf(this)
            % This function deterimines current user status is await, 
            % if current user already have rcf file in pool
            rcfpool = this.readall();
            tf = false;
            if isempty(rcfpool), return; end

            if this.current_await_only()
                tf = true;
                return;
            else
                for k = 1:numel(rcfpool)
                    if (rcfpool(k).status == "Await") ...
                        && (rcfpool(k).user_id == this.data.user_id)
                        tf = true;
                        return;
                    end
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
            elseif this.current_await_only()
                % only an old rcf exists and it is current user
                % set status to run, which will get resources in next 
                % "recycle" switch
                this.data.status = "Run";
            else
                uid = this.data.user_id;
                switch this.data.status
                    case "Await"
                        nprocgs = 0;                %number of processors that gpu shared
                        running_proc = false(numel(rcfpool), 1); % only acquired to running user

                        % generate task table:
                        % user  num_processor  t_estimate  resource  run_flag
                        tasks_ = table('Size', [numel(rcfpool), 5], ...
                            'VariableTypes', {'string','double','string','logical','logical'},...
                            'VariableNames', {'user', 'nproc', 'resrc', 'run', 'await'});

                        for k = 1:numel(rcfpool)
                            rcf_k = rcfpool(k);
                            tasks_.user(k) = rcf_k.user_id;
                            tasks_.nproc(k) = rcf_k.nworkers;
                            running_proc(k) = (rcf_k.progress > 0 && ...
                                rcf_k.progress < this.RECYCLE_REJECT_PROGRESS_THRESHOLD);
                            tasks_.resrc(k) = rcf_k.resource;
                            tasks_.run(k) = (rcf_k.status == "Run");
                            tasks_.await(k) = (rcf_k.status == "Await");
                        end

                        % which could cause a deadlock, todo
                        if all(~running_proc)
                            % until there exist a running user
                            % keep status: 'Await'
                            return;
                        end

                        % require workers
                        % elastic model:        cpu worker       gpu worker
                        %                [=====================|#############]
                        %                |<------------limit max------------>|
                        %                                      |<-limit max->|
                        req_user = [];
                        for k = 1:numel(rcfpool)
                            % require from the running user
                            if tasks_.run(k) == true
                                switch this.data.resource
                                    case "cpu"
                                        % cpu can only extension from cpu user
                                        if tasks_.resrc(k) == "cpu"
                                            req_user = [req_user; tasks_.user(k)]; %#ok<AGROW>
                                        else
                                            % protect gpu user, do not 
                                            % append them to required list
                                            nprocgs = nprocgs + tasks_.nproc(k);
                                        end
                                    case "cpu|gpu"
                                        % gpu user has higher extension
                                        % priority, no matter other users
                                        req_user = [req_user; tasks_.user(k)]; %#ok<AGROW>
                                    otherwise
                                end
                            end
                        end

                        % sum current user status for simple required
                        nruns = sum(tasks_.run);
                        if current_await_with_rcf(this)
                            nawaits = sum(tasks_.await);
                        else
                            % plus 1  for current user first join
                            nawaits = sum(tasks_.await) + 1;
                        end
                        
                        % return if too many user is running
                        if nruns >= this.USERS_NUMBER_MAX
                            this.await_counter = mod(this.await_counter+1, 5);
                            if this.await_counter == 1
                                waring("regrcf:tooManyUsersRunning", ...
                                    "Registration pool is busy, await for other tasks done...");
                            end

                            return;
                        else
                            this.await_counter = 0; % reset counter
                        end

                        % useful only at first requirement
                        % estimate average workers number
                        nproc_atot_req = max(nruns, ...
                            (this.nmax-nprocgs)/(nruns + nawaits));
                        
                        % required workers spread to running users
                        for k = 1:numel(req_user)
                            % even distribution resources
                            this.data.nworkers_req.(req_user(k)) ...
                                = [floor(nproc_atot_req / nruns), 0];
                        end

                        % check if the release sources are satisfied
                        nproc_atot_rel = 0;

                        for k = 1:numel(rcfpool)
                            if ismember(rcfpool(k).user_id, req_user)
                                % if required to current user
                                % extract the release user list
                                rel_user = string(fieldnames(rcfpool(k).nworkers_rel));
                                if ismember(uid, rel_user)
                                    % if current user in that release user list
                                    nproc_rel = rcfpool(k).nworkers_rel.(uid);
                                    
                                    % check if the resource was released
                                    % take the resource instead of current required
                                    % for avoiding infinity requirement
                                    if nproc_rel(2) == 1
                                        this.data.nworkers_req.(rcfpool(k).user_id) ...
                                            = nproc_rel;
                                        nproc_atot_rel = nproc_atot_rel + nproc_rel(1);
                                    end
                                end
                            end
                        end

                        % validate the sum of released resources and left 
                        % resource is enough
                        nlt = this.get_sys_nworkers_left();
                        if nlt + nproc_atot_rel >= floor(nproc_atot_req)
                            % resource is satisfied, update allocated workers, change status
                            this.data.nworkers = floor(nproc_atot_req);
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
                        rel_users = fieldnames(this.data.nworkers_rel);
                        % check others requirement
                        for k = 1:numel(rcfpool)
                            if rcfpool(k).status == "Await"
                                req_user = string(fieldnames(rcfpool(k).nworkers_req));
                                if ismember(uid, req_user)
                                    % if other user require to current user
                                    nproc_req = rcfpool(k).nworkers_req.(uid);

                                    % how many required, how many released
                                    % if this not mark released and that not mark
                                    % required
                                    if nproc_req(2) ~= 1 && ...
                                            ~ismember(rcfpool(k).user_id, rel_users)
                                        % set release flag to 1
                                        this.data.nworkers_rel.(rcfpool(k).user_id) ...
                                            = [nproc_req(1), 1];
                                        % let some people get rich first
                                        this.data.nworkers = ...
                                            this.data.nworkers - nproc_req(1);
                                        relsrc_ = relsrc_ + nproc_req(1);
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

            % if and only if all users are running and system has left
            % resoureces, allocate new workers
            switch this.data.status
                case "Run"
                    running_all = true;
                    nlt = this.get_sys_nworkers_left();
                    for k = 1:numel(rcfpool)
                        running_all = running_all & (rcfpool(k).status == "Run");
                    end

                    % average spread workers only if all users are running
                    dworkers = floor(nlt / numel(rcfpool));
                    if (dworkers > ceil(this.SPEED_UP_THRESHOLD*this.data.nworkers)) ...
                            && (running_all == true)
                        this.data.nworkers = this.data.nworkers + dworkers;
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

        function nlt = get_sys_nworkers_left(this)
            % This function reads pools and calculate left valid workers
            % note that n_max should be smaller than hardware supported maxima
            % for any user
            [~, n_max] = GetTaskWorkersMaxN(this.volopt, this.regopt);

            rcfpool = this.readall();

            if isempty(rcfpool)
                nlt = n_max;
            else
                nsp = 0;
                for k = 1:numel(rcfpool)
                    rcf_k = rcfpool(k);
                    if rcf_k.status == "Run"
                        nsp = nsp + rcf_k.nworkers;
                    end
                end
                nlt = n_max - nsp;
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