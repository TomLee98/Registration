classdef mQueue < matlab.mixin.Copyable
    %MQUEUE This class provide a basic queue and some operator
    
    properties (Access=private, Hidden)
        length;
    end

    properties(Access=private, Hidden, NonCopyable)
        data_v; % cell array is better choice
    end
    
    methods
        function this = mQueue()
            %MQUEUE Initialize an empty queue
            this.length = 0;
            this.data_v = {};
        end

        function enqueue(this, d)
            this.data_v = [this.data_v, {d}];
            this.length = this.length + 1;
        end

        function item = dequeue(this)
            if ~this.isempty()
                item = this.head();
                this.data_v(1) = [];
                this.length = this.length - 1;
            else
                item = [];
            end
        end

        function item = head(this)
            if ~this.isempty()
                item = this.data_v{1};
            else
                item = [];
            end
        end

        function item = tail(this)
            if ~this.isempty()
                item = this.data_v{end};
            else
                item = [];
            end
        end

        function tf = isempty(this)
            if this.size() == 0
                tf = true;
            else
                tf = false;
            end
        end

        function len = size(this)
            len = this.length;
        end

        function delete(this)
            this.data_v = {};
            this.length = 0;
        end
    end

    methods(Access = protected)
        function cpt = copyElement(this)
            cpt = copyElement@matlab.mixin.Copyable(this);
            cpt.data_v = this.data_v;   % cell copy
        end
    end
end

