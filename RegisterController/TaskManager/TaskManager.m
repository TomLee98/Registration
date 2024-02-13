classdef TaskManager < handle
    %TASKMANAGER This class is task manager object defination
    % 
    properties(Constant, Hidden)
        % ASK FOR REQUIRE DELAY
        REQUIRE_DELAY = 3;

        BUSY_RATIO = 2;
    end

    properties(Access=private, Hidden)
        regopts         % 1-by-1 regopt object, the registration options
        volopts         % 1-by-12 table, the volume series options
        regfrs          % 1-by-t positive integer, registration frame indices
        regedfn         % 1-by-1 positive integer, registered frames number
        movtmpl         % 1-by-1 regtmpl obejct, the registration template
        distrib         % 1-by-1 logical, distribution mode flag
        rcfobj          % 1-by-1 regrcf object, the content of resource communication file
                        % note that RCF should be a public protocal, shared
        taskqueue       % 1-by-1 mQueue, queue with Task object
        task_cur        % 1-by-1 Task object, the current task
        nworker_old     % 1-by-1 positive integer, the old usable workers number
        nworker_cur     % 1-by-1 positive integer, the current usable workers number
        parobj          % 1-by-1 parpool object

        sfolder         % 1-by-1 string, rcf pool folder
        nw_protected    % 1-by-1 double, the number of workers protected
        caller          % 1-by-1 Register obejct
    end

    properties(SetAccess=immutable, GetAccess=private, Hidden)
        PSMWN           % platform supported maximum workers number
                        % each user will get the same number
    end

    properties(GetAccess=?RegisterController,  Dependent)
        Task        % variable, get, auto update
    end

    properties(GetAccess=?Register, Dependent)
        SysInfo     % variable, get, manual update(intermittently require)
    end

    methods
        function this = TaskManager(volopt_, regopt_, movtmpl_, regfrs_)
            %TASKMANAGER A constructor
            arguments
                volopt_     (1,:)  table            % the movie information table
                regopt_     (1,1)  regopt           % the registration options struct
                movtmpl_    (1,1)  regtmpl          
                regfrs_     (1,:)  double {mustBePositive, mustBeInteger}
            end

            this.volopts = volopt_;
            this.regopts = regopt_;
            this.regfrs = regfrs_;
            this.movtmpl = movtmpl_;
            this.PSMWN = GetWorkersMaxN(this.regopts);
            this.nworker_old = 0;
            this.nworker_cur = 0;
            this.regedfn = 0;
        end

        function r = get.Task(this)
            r = this.task_cur;
        end

        function r = get.SysInfo(this)
            if this.distrib == true
                rcfpool = this.rcfobj.readall();

                r = table('Size', [numel(rcfpool), 6], ...
                    'VariableTypes', repmat({'string'},1,6), ...
                    'VariableNames',{'user','status','progress','time_used','n_cpu','n_gpu'});
                if numel(rcfpool) == 0
                    return;
                end

                % arrange by progress, sort as descend
                progress_ = nan(numel(rcfpool), 1);
                for n = 1:numel(rcfpool)
                    progress_(n) = str2double(rcfpool(n).progress(1));
                end
                [~, idx] = sort(progress_, "descend");
                for n = 1:numel(rcfpool)
                    m = idx(n);
                    time_used = string(datetime()-datetime(rcfpool(m).submit_time(1)));
                    r(n,:) = {rcfpool(m).user_id(1), rcfpool(m).status(1), rcfpool(m).progress(1), ...
                        time_used, rcfpool(m).nworkers(1), rcfpool(m).nworkers_max(1)};
                end
            else

            end
        end

        function setup(this, nwprotect_, distrib_, caller_, debug_)
            arguments
                this
                nwprotect_  (1,1)   double {mustBeNonnegative, mustBeInteger}
                distrib_    (1,1)   logical
                caller_     (1,1)   Register
                debug_      (1,1)   logical = false
            end
            this.nw_protected = nwprotect_;
            this.distrib = distrib_;
            this.caller = caller_;

            % link to rcf pool folder
            this.add_rcfpool();

            % communicate with rcf pool
            % here may hang out when no resource left
            % talk to will update workers number
            this.talk_to_rcfpool();
            
            % NOTE: update_taskqueue must running before update_parpool
            % update the task queue(depend on workers number)
            this.update_taskqueue();

            % update parpool
            % and adjust the current parpool workers number
            this.update_parpool(debug_);

            % get new task from taskqueue
            this.get_new_task();

            % set the main panel progress bar at 0
            this.caller.SetProgressBar(0);
        end

        function update(this, status_, debug_)
            % This function will update task queue:
            % if status is 0, set the task_cur status to "Done", enqueue 
            % task_cur (back to queue end)
            % if possible, dequeue front task of queue, if status is
            % "Await", update task_cur, else set task_cur to empty
            arguments
                this
                status_ (1,1)   double {mustBeInteger}
                debug_      (1,1)   logical = false
            end

            if status_ == 0
                % change the task running flag
                this.task_cur.Status = "Done";

                % push the task to the tail of taskqueue
                this.taskqueue.enqueue(this.task_cur);

                % update registration progress
                this.regedfn = this.regedfn + numel(this.task_cur.RegFrames);
                this.rcfobj.Progress = this.regedfn / numel(this.regfrs);
                this.caller.SetProgressBar(this.rcfobj.Progress);

                if this.distrib == true
                    % communicate with rcfs pool
                    this.talk_to_rcfpool();

                    % update the task queue(depend on workers number)
                    this.update_taskqueue();

                    % update parpool
                    this.update_parpool(debug_);
                end

                % get new task from taskqueue
                this.get_new_task();
            else
                throw(MException("TaskManager:update:invalidWorkerStatus", ...
                    "Can not create a new task. Workers status is %d.", status_));
            end
        end

        function delete(this)
            % free parpool resource
            delete(this.parobj);

            % remove rcf
            delete(this.rcfobj);
        end
    end

    methods(Access = private)

        function add_rcfpool(this)
            
            % setup rcf files folder
            if this.distrib == true
                if isunix()
                    sfpath = '/data/.rcfs';
                    success_flag = true;             % avoid users permission
                elseif ispc()
                    warning('off', 'MATLAB:MKDIR:DirectoryExists');
                    sfpath = 'C:\ProgramData\Reg3D\rcfs';
                    success_flag = mkdir(sfpath);    % generate folder except existed
                    warning('on', 'MATLAB:MKDIR:DirectoryExists');
                end
                if ~success_flag
                    throw(MException("TaskManager:mkdirFailed", ...
                        "Make resource communication folder failed."));
                end
                this.sfolder = string(sfpath);
                if ispc()
                    % hide the folder
                    fileattrib(this.sfolder, "+h");
                end
            else
                this.sfolder = "";  % no matter the sfolder, memory running
            end
            
        end

        function talk_to_rcfpool(this)
            % This is the most iomportant function in TaskManager, which
            % communicate with rcf pool, upload and download the schedule
            % update the current workers number
            if this.distrib == true
                if isempty(this.rcfobj) || ~isvalid(this.rcfobj)
                    % initialize an rcf object
                    this.rcfobj = regrcf(this.sfolder, this.volopts, this.regopts, this.regfrs);
                end

                % adjust if there exists rcf
                if this.rcfobj.is_pool_empty() ...
                        && (this.rcfobj.Status == "Await")
                    % if pool is empty and Status is "Await", 
                    % generate the first rcf, exclusive
                    this.rcfobj.NWorkersMax = max(this.PSMWN - this.nw_protected, 0);
                    this.rcfobj.NWorkers = this.rcfobj.NWorkersMax;

                    % require from rcf pool
                    this.rcfobj.update_resource("REQUIRE");
                    this.rcfobj.Status = "Run";
                    this.nworker_cur = ...
                            this.rcfobj.NWorkers - this.nw_protected;

                    % write to rcf pool
                    this.rcfobj.write();
                else
                    % if current user is await, wait for resource release
                    % omit memory or other hardwares limitation in
                    % distribution mode, TODO: more advanced
                    while this.rcfobj.Status == "Await"
                        % loop until rcf status changed to "Ready"
                        this.rcfobj.update_resource("REQUIRE");

                        this.rcfobj.write();

                        % delay some time
                        pause(this.REQUIRE_DELAY);
                    end
                    if this.rcfobj.Status == "Ready"
                        % new allocating branch
                        this.rcfobj.Status = "Run";
                        this.nworker_cur = ...
                            this.rcfobj.NWorkers - this.nw_protected;
                        this.rcfobj.write();
                    elseif this.rcfobj.Status == "Run"
                        % continue running branch
                        % any other is Await
                        if this.rcfobj.any_other_await()
                            % release some resource
                            this.rcfobj.update_resource("RELEASE");
                            this.nworker_cur = ...
                                this.rcfobj.NWorkers - this.nw_protected;

                            this.rcfobj.write();
                        else
                            % try to recycle some resources
                            this.rcfobj.update_resource("RECYCLE");
                            this.nworker_cur = ...
                                this.rcfobj.NWorkers - this.nw_protected;

                            this.rcfobj.write();
                        end
                    else
                        throw(MException("TaskManager:talk_to_rcfpool:innerError", ...
                            "Unexpected inner error: RCF status is invalid."));
                    end
                end
            else
                % exclusive mode
                % sfolder must be [], do nothing except update nworkers_cur
                this.nworker_cur = this.PSMWN - this.nw_protected;

                % update local rcf
                this.rcfobj = regrcf(this.sfolder, this.volopts, this.regopts, this.regfrs);
            end
        end

        function update_taskqueue(this)
            if this.nworker_cur ~= this.nworker_old
                % reset the task queue to match the workers
                % get all tasks which status are "Await", resort the frame
                % indices
                if ~isempty(this.taskqueue)
                    regframes = [];
                    while ~isempty(this.taskqueue)
                        task_ = this.taskqueue.dequeue();
                        if task_.Status == "Done"
                            break;
                        elseif task_.Status == "Await"
                            regframes = [regframes, task_.RegFrames]; %#ok<AGROW>
                        else
                            throw(MException("update_task:innerError", ...
                                "Invalid task status in task queue."));
                        end
                    end
                else
                    regframes = this.regfrs;
                end

                % use taskParser
                p = taskParser(this.movtmpl, regframes, this.volopts, ...
                    this.regopts, this.nworker_cur, this.distrib);
                p.parse(this.BUSY_RATIO);
                this.taskqueue = p.Results;
            end
        end

        function get_new_task(this)
            % get the next task in queue
            task_ = this.taskqueue.head();
            if ~isempty(task_) && (task_.Status == "Await")
                % new task dequeue
                this.task_cur = this.taskqueue.dequeue();
                % change the task running flag
                this.task_cur.Status = "Run";
            else
                this.task_cur = [];
            end
        end

        function update_parpool(this, debug_)
            if isempty(this.parobj)
                this.parobj = gcp("nocreate");
            end
            if ~isempty(this.nworker_cur) && ~isempty(this.nworker_old) ...
                    && (this.nworker_old~=this.nworker_cur)
                % close parpool if needed
                delete(this.parobj);

                % remove the possible latest failed jobs
                pcl = parcluster(parallel.defaultProfile);
                delete(pcl.Jobs);

                % make sure your cpu supports 'hyper-threads' technology
                % pcl.NumThreads = 2;

                % debug_ need spmdEnabled is true for mpiprofile and
                % mpiviewer working
                this.parobj = parpool(pcl, [1, this.nworker_cur], 'SpmdEnabled',debug_);

                % update workers number
                this.nworker_old = this.nworker_cur;
            end
        end
    end
end

function r = GetWorkersMaxN(regopt_)
switch regopt_.Mode
    case "global"
        r = GetCPUWorkersMaxN([], []);
    case "local"
        switch regopt_.Options.Hardware
            case "cpu"
                r = GetCPUWorkersMaxN([], []);
            case "cpu|gpu"
                r = GetGPUWorkersMaxN();
            otherwise
        end
    otherwise
end

end
