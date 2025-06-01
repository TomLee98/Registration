classdef clmhpr < handle
    %CLMHPR This class defines a clean helper to clean the mpimg/mpimgs
    % generated buffer file on local disk (on server as usual)

    properties(Access = private, Hidden)
        folder      (1,1)   string  % the clean folder
        strategy    (1,1)   struct  % the clean options
        cldaemon    (1,1)           % timer for activating cleaning 
    end

    properties(Access = public, Dependent)
        Capacity
        Trigger
    end
    
    methods
        function this = clmhpr(folder, capacity, trigger)
            %CLMHPR A Constructor
            arguments
                folder      (1,1)   string  {mustBeFolder}
                capacity    (1,1)   double  {mustBeInRange(capacity, 1, 1024)} = 128
                trigger     (1,1)   string  {mustBeMember(trigger, ["EXIT", "SIZE", "TIME", "OFF"])} = "EXIT"
            end

            this.folder = folder;

            this.strategy = struct("Capacity",  capacity, ...
                                   "Trigger",   trigger);

            this.cldaemon = timer("Period",60, "BusyMode","drop", ...
                "ExecutionMode","fixedRate","TasksToExecute",inf,...
                "Name","Timer_BufferCleanDaemon","TimerFcn",@this.CleanBuffer);
        end
        
        function r = get.Capacity(this)
            r = this.strategy.Capacity;
        end

        function set.Capacity(this, r)
            arguments
                this
                r   (1,1)   double  {mustBeInRange(r, 1, 1024)}
            end

            this.strategy.Capacity = r;
        end

        function r = get.Trigger(this)
            r = this.strategy.Trigger;
        end

        function set.Trigger(this, r)
            arguments
                this
                r   (1,1)   string  {mustBeMember(r, ["EXIT", "SIZE", "TIME", "OFF"])} = "EXIT"
            end

            this.strategy.Trigger = r;
        end

    end

    methods(Access = private)
        function CleanBuffer(this, ~, ~)
            % validate strategy and clean folder



        end
    end
end

