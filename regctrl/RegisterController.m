classdef RegisterController < handle
    %REGWORKER This class is regworker object defination
    % The RegWorker object can select registration algorithm and use the
    % TaskManager for multi-user calculation scheduling

    properties(Constant)
        KERNEL_OFF = 0
        KERNEL_ON = 1
    end

    properties(Access = private)
        regopts             % registration options
        nw_protected        % the protected workers number
        distrib             % distribution flag
    end

    properties(Access = private, Hidden)
        caller              % caller, must be Register object
        taskmgr             % TaskManager object
        regworker           % RegisterWorker object
        state               % 1-by-1 double, could be 0 or 1
        runtime             % 1-by-1 double, correction time used, seconds
    end

    properties(Access=public, Dependent)
        RegisterOptions     % variable, set,    running protected
        NWorkersProtected   % variable, set,    running protected
        Distributed         % variable, set,    running protected
    end

    properties(GetAccess=?Register, Dependent)
        State               % variable, get
        RunTime             % variable, get
        SysInfo             % variable, get
    end

    methods
        function this = RegisterController(caller_, regopt_, nwproctect_, distrib_)
            %REGISTER A constructor
            arguments
                caller_     (1,1)  Register
                regopt_     (1,1)  regopt
                nwproctect_ (1,2)  double {mustBeNonnegative, mustBeInteger} = [0, 0]
                distrib_    (1,1)  logical = false
            end

            this.caller = caller_;
            this.regopts = regopt_;
            this.nw_protected = nwproctect_;
            this.distrib = distrib_;

            % 
            this.taskmgr = TaskManager(regopt_);

            this.state = this.KERNEL_OFF;
        end

        function set.RegisterOptions(this, r_)
            arguments
                this
                r_  (1,1)   regopt
            end

            if this.state == this.KERNEL_OFF
                this.regopts = r_;
            else
                warning("RegisterController:invalidOperation", ...
                    "Can not set registration options when kernel is running.");
            end
        end

        function set.NWorkersProtected(this, r_)
            arguments
                this
                r_  (1,2)   double  {mustBeNonnegative, mustBeInteger}
            end

            if this.state == this.KERNEL_OFF
                this.nw_protected = r_;
            else
                warning("RegisterController:invalidOperation", ...
                    "Can not set workers number when kernel is running.");
            end
        end

        function r = get.NWorkersProtected(this)
            r = this.nw_protected;
        end

        function set.Distributed(this, r_)
            arguments
                this
                r_  (1,1)   logical
            end

            if this.state == this.KERNEL_OFF
                this.distrib = r_;
            else
                warning("RegisterController:invalidOperation", ...
                    "Can not set distributed state when kernel is running.");
            end
        end

        function r = get.State(this)
            r = this.state;
        end

        function r = get.RunTime(this)
            r = round(this.runtime, 1);
        end

        function r = get.SysInfo(this)
            r = struct("T",     this.taskmgr.SysInfo, ...
                       "NWP",   this.nw_protected);
        end

    end

    methods(Access=public)
        function run(this, movraw_, movaligned_, movtmpl_, regfr_)
            % This function is the controller of registration
            arguments
                this
                movraw_     (1,1)   regmov
                movaligned_ (1,1)   regmov
                movtmpl_    (1,1)   regtmpl
                regfr_      (1,:)   double {mustBePositive, mustBeInteger}
            end

            % setup the task manager
            this.taskmgr.setup(this.regopts, movraw_, movtmpl_, regfr_, ...
                this.nw_protected, this.distrib, this.caller);

            % initialize register worker
            this.regworker = RegisterWorker(movraw_, movaligned_);

            % turn on engine
            this.state = this.KERNEL_ON;
            st = tic;

            % loop until no valid task
            task = this.taskmgr.Task;
            while ~isempty(task)
                if this.state == this.KERNEL_OFF
                    break;
                end

                status = this.regworker.correct(task);

                this.taskmgr.update(status);    % could hang out for resources
                
                task = this.taskmgr.Task;
            end

            this.runtime = toc(st);

            this.clear_after_running(movraw_, regfr_);
        end

        function stop(this)
            % change state
            this.state = this.KERNEL_OFF;
        end

        function delete(this)
            delete(this.taskmgr);
            delete(this.regworker);

            % ~
        end
    end

    methods(Access=private, Hidden)
        function clear_after_running(this, movraw_, regfr_)
            % clean resource
            clear(this.taskmgr);
            delete(this.regworker);

            switch this.state
                case this.KERNEL_ON      % normal exit
                    % send the complete message
                    switch this.regopts.Algorithm
                        case "MANREG"
                            this.caller.SendMessage(["REGISTRATION_FRAMES", ...
                                sprintf(" %.1f", movraw_.Time(regfr_))]);
                        otherwise
                            this.caller.SendMessage(["REGISTRATION_COMPLETE",""]);
                            this.caller.SendMessage(["TIME_COST_SEC", sprintf(" %.1f", this.RunTime)]);
                    end

                    this.state = this.KERNEL_OFF;
                case this.KERNEL_OFF     % force exit
                    % send the stop message
                    this.caller.SendMessage(["TERMINATE_REGISTRATION",""]);
                    this.caller.SendMessage(["TIME_COST_SEC", sprintf(" %.1f", this.RunTime)]);
                otherwise
            end
        end
    end
end
