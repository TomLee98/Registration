classdef OperationHistoryManager < handle
    %REGHIMGR This class defines OperationHistoryManager for reg3D app, 
    % which is the agent to handle image and operations resources
    % NOTE: change storage strategy need

    properties(Access = public, Dependent)
        Strategy        % ___/get, 1-by-1 string, could be "PERFORMANCE"/"RESOURCE"/"BALANCE"
        MemoryUsed      % ___/get, 1-by-1 double, unit as GBytes
        DiskUsed        % ___/get, 1-by-1 double, unit as GBytes
        CurrentData     % ___/get, 1-by-1 regmov, indicating current activated node data
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
        branch_cur                          % current operated branch
        node_active      (1,1)              % current active node
    end
    
    methods
        function r = get.CurrentData(this)
            if this.isempty()
                r = regmov.empty();
            else
                r = getData(this.node_active.RSPoint);
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
        function this = OperationHistoryManager(op_tree, op_txt, ctxt_menu, strategy)
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
            end

            this.optree = op_tree;
            this.optree.SelectionChangedFcn = @(~, event)this.sltchg_callback(event);
            this.optxt = op_txt;
            this.ctmenu = ctxt_menu;
            this.strategy = strategy;

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
                    node_data.Properties.Location = file_dp.File;
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
            ActivateNode(this, node_new);
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

            if isequal(this.node_active, this.optree.SelectedNodes)
                disp("Node is already activated!");
            else
               %% Deactivated current active node
               this.deactivate();

               %% Activate selected node
               this.activate(node);
            end
            
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

        % This is select changed callback binding on tree
        function sltchg_callback(this, event)
            nd = event.SelectedNodes.NodeData;

            if ~isempty(nd)
                %% parse node data and generate formatted text
                txt = sprintf("[Time] %s\n[Operation] %s\n[Properties]{\n%%s}\n", ...
                    string(nd.Time), nd.Operation);
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
                                    end
                                case "Frames"
                                    txt_prop = txt_prop + sprintf("\t[%s] %s\n", string(props{k}), ...
                                            string(nd.Properties.(props{k})).join(":").join(","));
                                otherwise
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
end

