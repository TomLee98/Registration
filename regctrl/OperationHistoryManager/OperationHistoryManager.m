classdef OperationHistoryManager < handle
    %REGHIMGR This class defines OperationHistoryManager for reg3D app, 
    % which is the agent to handle image and operations resources
    % NOTE: change storage strategy need

    properties(Access = public, Dependent)
        Strategy        % ___/get, 1-by-1 string, could be "PERFORMANCE"/"RESOURCE"/"BALANCE"
        IsDistributed   % set/get, 1-by-1 logical, indicates if data spread on disk
        MemoryUsed      % ___/get, 1-by-1 double, unit as GBytes
        DiskUsed        % ___/get, 1-by-1 double, unit as GBytes
        CurrentData     % ___/get, 1-by-1 regrspt, indicating current activated restored node
        CurrentNodeTag  % ___/get, 1-by-1 string, could be empty
    end

    properties(SetAccess = immutable, GetAccess = private)
        optree
        optxt
        ctmenu
        strategy
        flagsrc = string(fileparts(mfilename("fullpath"))).extractBefore(...
            "regctrl") + filesep + "sources" + filesep + "active.png"
        readysrc = string(fileparts(mfilename("fullpath"))).extractBefore(...
            "regctrl") + filesep + "sources" + filesep + "ready.png"
    end

    properties(Access = private, Hidden)
        branch_cur                              % current operated branch
        node_active     (1,1)                   % current active node
        is_distrib      (1,1)   logical = false % data distributed flag 
        n_nodes         (1,1)   double  = 0     % number of nodes in a tree
    end
    
    methods
        function r = get.CurrentData(this)
            if this.isempty()
                r = [];
            else
                r = this.node_active.RSPoint;
            end
        end

        function r = get.CurrentNodeTag(this)
            r = string(this.node_active.Tag);
        end

        function r = get.IsDistributed(this)
            r = this.is_distrib;
        end

        function set.IsDistributed(this, r)
            arguments
                this 
                r       (1,1)   logical
            end

            % try to spread or gather all data
            if r == true && this.is_distrib == false
                % spread data to disk
                this.move_storage("spread");

                this.is_distrib = true;
            elseif r == false && this.is_distrib == true
                % gather data to memory
                this.move_storage("gather");

                this.is_distrib = false;
            end
        end

        function r = get.Strategy(this)
            r = this.strategy;
        end

        function r = get.MemoryUsed(this)
            r = nan;
        end

        function r = get.DiskUsed(this)
            r = nan;
        end
        
    end

    methods (Access = public)
        function this = OperationHistoryManager(op_tree, op_txt, ctxt_menu, strategy, dbflag)
            %REGHIMGR A Constructor
            % Input:
            %   - op_tree: 1-by-1 matlab.ui.container.Tree as nodes container
            %   - op_txt: 1-by-1 matlab.ui.container.TextArea as
            %               information viewer
            %   - ctxt_menu: 1-by-1 matlab.ui.container.ContextMenu as
            %               interaction callback adapter
            %   - strategy: 1-by-1 string, could be "performance"/"resource"/"balance"
            arguments
                op_tree     (1,1)   matlab.ui.container.Tree
                op_txt      (1,1)   matlab.ui.control.TextArea
                ctxt_menu   (1,1)   matlab.ui.container.ContextMenu
                strategy    (1,1)   string  {mustBeMember(strategy, ...
                                    ["PERFORMANCE", "RESOURCE", "BALANCE"])} = "PERFORMANCE"
                dbflag      (1,1)   logical = false
            end

            this.optree = op_tree;
            this.optree.SelectionChangedFcn = @(~, event)this.sltchg_callback(event);
            this.optxt = op_txt;
            this.ctmenu = ctxt_menu;
            this.strategy = strategy;
            this.is_distrib = dbflag;

            this.branch_cur = mStack();

            this.node_active = this.optree;
        end

        function delete(this)
            % free all data


        end

        function tf = isempty(this)
            % only tree without any node
            tf = isa(this.node_active, "matlab.ui.container.Tree");
        end

        function AddNode(this, optr, args, varargin)
            % This function create new node at optree and store a regrspt
            % object in the node
            %% create restore point object
            switch this.Strategy
                case "PERFORMANCE"
                    rs_node = regrspt(optr, args, varargin{:});
                case "RESOURCE"
                    warning("OperationHistoryManager:unfinishedItem", ...
                        "Developing...");
                case "BALANCE"
                    warning("OperationHistoryManager:unfinishedItem", ...
                        "Developing...");
                otherwise
                    throw(MException("OperationHistoryManager:invalidStrategy", ...
                        "Unsupported transmission optimization strategy."));
            end

            %% create node on the op tree
            node_data = struct("Operation",  [], ...
                               "Properties", [], ...
                               "RSPoint",    rs_node, ...
                               "Time",       datetime("now", "Format","MM/dd HH:mm:ss"));
            switch optr
                case constdef.OP_LOAD
                    node_text = "load";
                    file_dp = rs_node.File;
                    % node_data.Properties.Location = file_dp.File;
                    node_data.Properties.Dims = rs_node.ImageDim;
                    switch file_dp.File
                        case "Memory"
                            node_data.Properties.Size = round(file_dp.Size.mem/2^30, 1);
                        otherwise
                            node_data.Properties.Size = round(file_dp.Size.map/2^30, 1);
                    end
                    
                case constdef.OP_CROP
                    node_text = sprintf("crop(%s)", rs_node.CropDim);

                    % append individual field
                    switch rs_node.CropDim
                        case "XY"
                            node_data.Properties = struct("Origin", args.xy(1, 1:2), ...
                                                          "Width",  args.xy(2, 1), ...
                                                          "Height", args.xy(2, 2));
                        case "Z"
                            node_data.Properties = struct("Origin", args.z(1), ...
                                                          "Slices", args.z(2)-args.z(1)+1);
                        case "T"
                            frs = calc_blocked_t(args);
                            node_data.Properties = struct("Frames", frs);
                        otherwise
                    end

                case constdef.OP_REGISTER
                    node_text = sprintf("register(%s)", args.Arguments.Mode);
                    node_data.Properties = struct("Algorithm",      args.Arguments.Algorithm, ...
                                                  "SubAlgorithm",   args.Arguments.SubAlgorithm);

                case constdef.OP_SEGMENT
                    node_text = sprintf("segment(%d)", rs_node.CellsCount);
                    node_data.Properties = struct("Method", args.method);

                otherwise
            end
            % append operation field
            node_data.Operation = optr;

            %% activate new node
            % change the current active node
            % modify: Text, Icon, NodeData, ContextMenu
            node_new = uitreenode(this.node_active, "Text",node_text, ...
                "NodeData",node_data, "ContextMenu",this.ctmenu);

            if ~this.isempty(), this.node_active.expand(); end % expand

            ActivateNode(this, node_new);                       % activate new node

            this.n_nodes = this.n_nodes + 1;
            this.node_active.Tag = string(this.n_nodes);

            %% 

            function r = calc_blocked_t(arg)
                frames = arg.f;
                r = zeros(0, 2);
                if isscalar(frames)
                    r(1) = frames; r(2) = frames;
                    return;
                end

                % multi frames selected
                ps = 1;
                pe = 2;
                while pe <= numel(frames)
                    if pe < numel(frames)
                        if frames(pe) - frames(pe-1) == 1
                            pe = pe + 1;
                        else
                            r = [r; frames(ps), frames(pe-1)]; %#ok<AGROW>
                            ps = pe;
                            pe = pe + 1;
                        end
                    else    % EOF
                        r = [r; frames(ps), frames(pe)]; %#ok<AGROW>
                        break;
                    end
                end
            end
        end

        function DelNode(this, node)
            arguments
                this 
                node (1,1)  matlab.ui.container.TreeNode
            end

            % 
            disp("Delete Test!");
        end

        function ActivateNode(this, node)
            arguments
                this 
                node (1,1)  matlab.ui.container.TreeNode
            end

            if isequal(this.node_active, node)
                disp("Node is already activated!");
            else
               % Deactivated current active node
               this.deactivate();

               % select the node
               this.optree.SelectedNodes = node;

               % activate selected node
               this.activate(node);

               % refresh snapshot appearance
               this.refresh_snapshot();
            end
            
        end

        % This function selects next node after current selection
        function SelectNextNode(this)
            % find the node with tag is the next tag

        end

        % This function gets one previous node before given node, if given
        % node is empty, active node replaced
        function node = GetPreviousNodeBefore(this, node)
            arguments
                this 
                node (1,:)  matlab.ui.container.TreeNode = []
            end

            if isempty(node)
                % return node before active node
                if isa(this.node_active.Parent, "matlab.ui.container.TreeNode")
                    node = this.node_active.Parent;
                else
                    node = [];
                end
            else
                if isa(node.Parent, "matlab.ui.container.TreeNode")
                    node =  node.Parent;
                else
                    node = [];
                end
            end
        end

        % This function gets all operations at given node
        function optrs = GetAllPreviousOperatorsAt(this, node)
            arguments
                this 
                node (1,:)  matlab.ui.container.TreeNode = []
            end

            
        end

        function PackageProjectTo(this, folder)
            % This function package all data as project to given folder

        end

        function LoadProjectFrom(this, folder)
            % This function load data from given project folder

        end

    end

    methods (Access = private)

        % This function activates node
        function activate(this, node)
            % node: 1-by-1 matlab.ui.container.TreeNode

            this.node_active = node;
            this.node_active.Icon = this.flagsrc;
        end


        % This function deactivates current active node
        function deactivate(this)
            % deactivate node in nonempty tree
            if ~this.isempty()
                % deactivate means
                % [1] drop (all )branch related to active node accord to strategy
                % [2] set the active node to tree (root)
                % [3] update some flag

                % remove the older node active flag
                this.node_active.Icon = this.readysrc;

                this.node_active = this.optree;
            end

        end

        % This function moves the data storage between memory and disk
        function move_storage(this, optr)
            arguments
                this
                optr    (1,1)   string  {mustBeMember(optr, ["spread", "gather"])}
            end

            % use depth-first traversal to spread data onto disk
            st = mStack();
            st.push(this.optree);   % push the tree root
            while ~isempty(st)
                % pop the node
                node = st.pop();

                % visit node
                if isa(node, "matlab.ui.container.TreeNode")
                    nd = node.NodeData;
                    switch optr
                        case "spread"
                            nd.RSPoint.IsDistributed = true;    % spread data
                        case "gather"
                            nd.RSPoint.IsDistributed = false;    % spread data
                        otherwise
                    end
                end

                nodes = node.Children;     % nodes
                for k = 1:numel(nodes), st.push(nodes(k)); end
            end
        end

        % This is select changed callback binding on tree
        function sltchg_callback(this, event)
            % get node data
            node = event.SelectedNodes;

            % refresh snapshot appearance
            this.refresh_snapshot(node);
        end

        function refresh_snapshot(this, node)
            arguments
                this
                node    (1,:)   = []
            end

            % set default node as active node
            if isempty(node), node = this.node_active; end
            nd = node.NodeData;

            %% parse node data and generate formatted text
            txt = sprintf("[Tag] %s\n[Time] %s\n[Operation] %s\n[Properties]{\n%%s}\n", ...
                node.Tag, string(nd.Time), nd.Operation);
            txt_prop = "";
            props = fieldnames(nd.Properties);
            switch nd.Operation
                case constdef.OP_LOAD
                    for k = 1:numel(props)
                        % key-value as props{k}-nd.Properties.(props{k})
                        switch props{k}
                            case "Dims"
                                txt_prop = txt_prop + sprintf("\t[%s] (%s)\n", string(props{k}), ...
                                    string(nd.Properties.(props{k})).join(","));
                            case "Size"
                                txt_prop = txt_prop + sprintf("\t[%s] %s GBytes\n", string(props{k}), ...
                                    string(nd.Properties.(props{k})));
                            otherwise
                                txt_prop = txt_prop + sprintf("\t[%s] %s\n", string(props{k}), ...
                                    string(nd.Properties.(props{k})));
                        end
                    end
                case {constdef.OP_REGISTER, constdef.OP_SEGMENT}
                    for k = 1:numel(props)
                        % key-value as props{k}-nd.Properties.(props{k})
                        txt_prop = txt_prop + sprintf("\t[%s] %s\n", string(props{k}), ...
                            string(nd.Properties.(props{k})));
                    end
                case constdef.OP_CROP
                    for k = 1:numel(props)
                        % key-value as props{k}-nd.Properties.(props{k})
                        switch props{k}
                            case "Origin"
                                if isequal("XY", nd.RSPoint.CropDim)
                                    txt_prop = txt_prop + sprintf("\t[%s] (%s)\n", string(props{k}), ...
                                        string(nd.Properties.(props{k})).join(","));
                                else
                                    txt_prop = txt_prop + sprintf("\t[%s] %s\n", string(props{k}), ...
                                        string(nd.Properties.(props{k})));
                                end
                            case "Frames"
                                txt_prop = txt_prop + sprintf("\t[%s] %s\n", string(props{k}), ...
                                    string(nd.Properties.(props{k})).join(":").join(","));
                            otherwise
                                txt_prop = txt_prop + sprintf("\t[%s] %s\n", string(props{k}), ...
                                    string(nd.Properties.(props{k})));
                        end
                    end
                otherwise
            end

            txt = sprintf(txt, txt_prop);   % nested generation

            %% update text field to display the snapshot
            this.optxt.Value = txt;
        end
    end
end

