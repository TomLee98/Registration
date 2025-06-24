classdef regohm < handle
    %REGOHM This class defines "Operation History Manager" for reg3D app, 
    % which is the agent to handle image and operations resources

    properties(Constant, Hidden)
        OP_SKIP_OUTPUT = {constdef.OP_LOAD, constdef.OP_SEGMENT}
    end

    properties(Access = public, Dependent)
        ActiveNodeData  % ___/get, 1-by-1 struct with field {arg, cbot, data, optr}
        ActiveNodeTag   % ___/get, 1-by-1 string, could be empty
        IsDistributed   % set/get, 1-by-1 logical, indicates if data spread on disk
        CachePolicy     % ___/get, 1-by-1 string, could be "PERFORMANCE"/"RESOURCE"/"BALANCE"
    end

    properties(SetAccess = immutable, GetAccess = private)
        ufig            % 1-by-1 matlab.ui.Figure
        optree          % 1-by-1 matlab.ui.container.Tree
        optxt           % 1-by-1 matlab.ui.controller.TextArea
        ctmenu          % 1-by-1 matlab.ui.container.ContextMenu
        optsty          % 1-by-1 string, indicate the manager running method
        flagsrc = string(fileparts(mfilename("fullpath"))).extractBefore(...
            "regctrl") + filesep + "sources" + filesep + "active.png"
        readysrc = string(fileparts(mfilename("fullpath"))).extractBefore(...
            "regctrl") + filesep + "sources" + filesep + "ready.png"
        calcsrc = string(fileparts(mfilename("fullpath"))).extractBefore(...
            "regctrl") + filesep + "sources" + filesep + "rtcalc.png"
    end

    properties(Access = private, Hidden)
        dptr            (1,1)   regmov = regmov.empty()     % current node related data pointer
        is_distributed  (1,1)   logical = false             % data distributed flag 
        node_active     (1,1)                               % current active node (or tree)
        nodes_total_num (1,1)   double  = 0                 % number of total generated nodes in the tree
    end
    
    methods
        function r = get.ActiveNodeData(this)
            if this.isempty()
                r = struct("arg",   [], ...
                           "cbot",  [], ...
                           "data",  this.dptr, ...
                           "others",[], ...
                           "optr",  "");
                return;
            end

            % load data in the current restore point
            rsp = this.node_active.NodeData.RSPoint;

            r = struct("arg",   rsp.Arguments, ...
                       "cbot",  rsp.Segmentor, ...
                       "data",  this.dptr, ...
                       "others",rsp.Others, ...
                       "optr",  rsp.Operator);
        end

        function r = get.ActiveNodeTag(this)
            r = string(this.node_active.Tag);
        end

        function r = get.IsDistributed(this)
            r = this.is_distributed;
        end

        function set.IsDistributed(this, r)
            arguments
                this 
                r       (1,1)   logical
            end

            % try to spread or gather all data
            if r == true && this.is_distributed == false
                % spread data to disk
                this.move_storage("spread");

                this.is_distributed = true;
            elseif r == false && this.is_distributed == true
                % gather data to memory
                this.move_storage("gather");

                this.is_distributed = false;
            end
        end

        function r = get.CachePolicy(this)
            r = this.optsty;
        end
        
    end

    methods (Access = public)
        function this = regohm(ufig, op_tree, op_txt, ctxt_menu, strategy, distributed)
            %REGHIMGR A Constructor
            % Input:
            %   - ufig: 1-by-1 matlab.ui.Figure for progress bar display
            %   - op_tree: 1-by-1 matlab.ui.container.Tree as nodes container
            %   - op_txt: 1-by-1 matlab.ui.container.TextArea as
            %               information viewer
            %   - ctxt_menu: 1-by-1 matlab.ui.container.ContextMenu as
            %               interaction callback adapter
            %   - strategy: 1-by-1 string, could be "performance"/"resource"/"balance"
            %   - distributed: 1-by-1 logical, indicate if is distributed mode
            arguments
                ufig        (1,1)   matlab.ui.Figure
                op_tree     (1,1)   matlab.ui.container.Tree
                op_txt      (1,1)   matlab.ui.control.TextArea
                ctxt_menu   (1,1)   matlab.ui.container.ContextMenu
                strategy    (1,1)   string  {mustBeMember(strategy, ...
                                    ["PERFORMANCE", "RESOURCE", "BALANCE"])} = "PERFORMANCE"
                distributed (1,1)   logical = false
            end

            this.ufig = ufig;
            this.optree = op_tree;
            this.optree.SelectionChangedFcn = @(~, event)this.sltchg_callback(event);
            this.optxt = op_txt;
            this.ctmenu = ctxt_menu;
            this.optsty = strategy;
            this.is_distributed = distributed;
            this.node_active = this.optree;
        end

        function delete(this)
            % free all data
            % from the newest node to oldest node

            if ~isempty(this.optree) && isvalid(this.optree) ...
                    && ~isempty(this.optree.Children)
                nodes = this.optree.Children;
                for k = 1:numel(nodes)
                    this.remove(nodes(k));
                end
            end
        end

        function tf = isempty(this)
            % only tree without any node
            tf = isa(this.node_active, "matlab.ui.container.Tree");
        end

        function AddNode(this, optr, args, varargin)
            % This function create new node at optree and store a regrspt
            % object in the node

            %% append node after active node
            node_new = this.create_new_node(optr, args, varargin);

            %% activate new node
            this.ActivateNode(node_new);
        end

        function DelNode(this, node)
            arguments
                this 
                node (1,1)  matlab.ui.container.TreeNode
            end

            %% auto change the active node
            this.automove_active_node(node);

            %% free memory/disk
            this.remove(node);

            %% update appearance
            this.update_manage_view();
        end

        function ActivateNode(this, node)
            arguments
                this 
                node (1,1)  matlab.ui.container.TreeNode
            end

            if ~isequal(this.node_active, node)
               %% deactivated current active node
               this.deactivate();

               %% activate given node
               this.optree.SelectedNodes = node;
               this.activate(node);

               %% update appearance
               this.update_manage_view();
            end
        end

        % This function gets all operations at given node
        function optrs = GetAllPreviousOperatorsAt(this, node)
            arguments
                this 
                node (1,:) = []
            end

            if this.isempty(), optrs = string([]); return; end

            if isempty(node), node = this.node_active; end

            optrs = strings([]);

            % use stack for previous tracing
            nodes = mStack();
            nodes.push(node);   % push the first node
            while ~isempty(this.get_node_before(node))
                node = this.get_node_before(node);
                nodes.push(node);
            end

            % pop node should be "append order"
            while ~isempty(nodes)
                node = nodes.pop();
                switch node.NodeData.RSPoint.Operator
                    case this.OP_SKIP_OUTPUT
                        % omit the operation
                    case constdef.OP_CROP
                        switch node.NodeData.RSPoint.CropDim
                            case {"XY", 'Z'}
                                optrs = [optrs, "crop"]; %#ok<AGROW>
                            case "T"
                                optrs = [optrs, "cut"]; %#ok<AGROW>
                            otherwise
                        end
                    case constdef.OP_REGISTER
                        optrs = [optrs, "aligned"]; %#ok<AGROW>
                    otherwise
                end
            end
        end

        % This function selects next node after current selection
        function SelectNextNode(this, loop)
            arguments
                this 
                loop    (1,1)   logical  = false
            end

            %% take next node
            node = this.optree.SelectedNodes;
            node = this.get_node_after(node);
            if ~isempty(node)
                this.optree.SelectedNodes = node;
            else
                if loop == true
                    % start at first
                    this.optree.SelectedNodes = this.optree.Children(1);
                else
                    % skip, hold on
                end
            end

            %% update appearance
            this.update_snapshot_view();    % snap only for fast view
        end

        % This function packages all data as project to given folder
        function PackageProjectTo(this, folder)
            

        end

        % This function loads data from given project folder
        function LoadProjectFrom(this, folder)
            

        end

    end

    methods (Access = private)

        % This function creates a new node as children of active node
        function node = create_new_node(this, optr, args, vars)

            %% create restore point object by cache policy
            is_storage = constdef.TRANSMISSION_STORAGE_TABLE{this.optsty, optr};

            if is_storage == false
                % drop data node and replace by empty
                for n = 1:2:numel(vars)
                    if isequal("Data", vars{n})
                        vars{n+1} = regmov.empty();     % drop, handle only be holded by Register
                        break;
                    end
                end
            end

            % construct restore point
            rs_node = regrspt(optr, args, vars{:});

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
                    st = args.(constdef.OP_CROP);

                    % append individual field
                    switch rs_node.CropDim
                        case "XY"
                            node_data.Properties = struct("Origin", st.xy(1, 1:2), ...
                                                          "Width",  st.xy(2,1)-st.xy(1,1)+1, ...
                                                          "Height", st.xy(2,2)-st.xy(1,2)+1);
                        case "Z"
                            node_data.Properties = struct("Origin", st.z(1), ...
                                                          "Slices", st.z(2)-st.z(1)+1);
                        case "T"
                            frs = regohm.reformat_time(st);
                            node_data.Properties = struct("Frames", frs);
                        otherwise
                    end
                case constdef.OP_REGISTER
                    st = args.(constdef.OP_REGISTER);
                    %TODO: add more descriptor
                    node_text = sprintf("register(%s)", st.Mode);
                    switch st.Algorithm
                        case {"TCREG", "OCREG"}
                            switch st.Mode
                                case "global"
                                    node_data.Properties = struct("Algorithm",      st.Algorithm, ...
                                                                  "Transformation", st.Options.TformType, ...
                                                                  "CoarseRegister", st.Options.CoarseAlg, ...
                                                                  "MaxStep",        string(st.Options.MaxStep), ...
                                                                  "MinStep",        string(st.Options.MinStep));
                                case "local"
                                    node_data.Properties = struct("Algorithm",          st.SubAlgorithm, ...
                                                                  "Compensation",       string(st.Options.ImageRehist), ...
                                                                  "CompensationGrade",  string(st.Options.RepAcc));
                                otherwise
                            end
                        case "MANREG"
                            node_data.properties = struct("Transformation", st.Options.TformType, ...
                                                          "Degree",         string(st.Options.Degree), ...
                                                          "ProjectedView",  st.Options.DView, ...
                                                          "Isometric",      string(st.Options.Isometric));
                        case "LTREG"
                            node_data.Properties = struct("Keyframes",      string(reshape(st.Options.Keyframes,1,[])).join(","), ...
                                                          "AutoKeyframes",  string(st.Options.AutoKeyframe), ...
                                                          "MaxStep",        string(st.Options.MaxStep), ...
                                                          "MinStep",        string(st.Options.MaxStep));
                        otherwise
                    end
                case constdef.OP_SEGMENT
                    st = args.(constdef.OP_SEGMENT);
                    node_text = sprintf("segment(%d)", rs_node.CellsCount);
                    node_data.Properties = struct("Method", st.method);
                otherwise
            end

            node_data.Operation = optr;

            %% add node to active node
            node = uitreenode(this.node_active, "Text",node_text, ...
                "NodeData",node_data, "ContextMenu",this.ctmenu);
            this.nodes_total_num = this.nodes_total_num + 1;
            node.Tag = string(this.nodes_total_num);
        end

        % This function activates node
        function activate(this, node)
            % node: 1-by-1 matlab.ui.container.TreeNode

            % activate node (get data) by strategy
            rpt = node.NodeData.RSPoint;

            is_storage = constdef.TRANSMISSION_STORAGE_TABLE{this.optsty, rpt.Operator};

            if is_storage == true
                % just read from node
                this.dptr = rpt.getData();
            else
                % find the shortest calculation path, equals to find the 
                % nearest storage which is of a higher  seniority in the 
                % family hierarchy

                % use stack for previous tracing
                nodes = mStack();
                nodes.push(node);   % push the first node
                while ~isempty(this.get_node_before(node))
                    % push node
                    node = this.get_node_before(node);
                    nodes.push(node);

                    % determine breaking if node strategy is storage
                    rpt = node.NodeData.RSPoint;
                    is_storage = constdef.TRANSMISSION_STORAGE_TABLE{this.optsty, rpt.Operator};
                    if is_storage, break; end
                end

                % generate calculation waitbar
                msg = sprintf("Calculating ...(0/%d)", nodes.size()-1);
                udlg = uiprogressdlg(this.ufig, "Indeterminate","off", ...
                    "Title","Calculation Progress", "Message",msg, ...
                    "Value",0, "Icon", this.calcsrc);

                % real time calculation chain
                node = nodes.pop();
                rpt = node.NodeData.RSPoint;
                this.dptr = rpt.getData();     % get storage
                n = 0; n0 = nodes.size();
                while ~isempty(nodes)
                    node = nodes.pop();
                    rpt = node.NodeData.RSPoint;

                    % update current dptr
                    this.dptr = rpt.getData(this.dptr);

                    n = n + 1;
                    udlg.Value = n / n0;
                    udlg.Message = sprintf("Calculating ...(%d/%d)", n, n0);
                    drawnow nocallbacks;
                end

                delete(udlg);
            end

            % update active node
            this.node_active = node;
        end

        % This function deactivates current active node
        function deactivate(this)
            % deactivate node in nonempty tree
            if ~this.isempty()
                if isMATLABReleaseOlderThan("R2024b")
                    % modify the flag
                    this.node_active.Icon = '';
                end
                % change active node as tree
                this.node_active = this.optree;
            end

        end

        % This function auto moves the active node if it is on the branch
        % which will be deleted
        function automove_active_node(this, node)
            if this.is_posterity_of(node)
                % priority: brother > parent
                if numel(node.Parent.Children) > 1
                    brothers = setdiff(node.Parent.Children, node);
                    this.ActivateNode(brothers(1));
                else
                    if isa(node.Parent, "matlab.ui.container.TreeNode")
                        % activate the nearest
                        this.ActivateNode(node.Parent);
                    else
                        % skip, tree only
                    end
                end
            end
        end

        % This function update manager appearance accords to current 
        % tree status
        function update_manage_view(this)
            % global control
            this.update_operation_view();
            this.update_source_view();

            % local control
            this.update_snapshot_view();

            drawnow
        end

        function remove(this, node)
            arguments
                this %#ok<INUSA>
                node (1,1)  matlab.ui.container.TreeNode
            end

            nodes_delete = mStack();
            nodes_visit = mStack();

            nodes_visit.push(node);

            while ~isempty(nodes_visit)
                % push to stack
                node = nodes_visit.pop();
                nodes_delete.push(node);

                % en-queue children
                nodes = node.Children;
                for k = 1:numel(nodes), nodes_visit.push(nodes(k));end
            end

            % delete on each node
            while ~isempty(nodes_delete)
                node = nodes_delete.pop();

                % invoke delete
                delete(node.NodeData.RSPoint);
            end

            % remove node handle
            delete(node);
        end

        % This function gets one previous node before given node, if given
        % node is empty, active node replaced
        function node = get_node_before(this, node)
            arguments
                this 
                node (1,:) = []
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

        function node = get_node_after(this, node, policy)
            arguments
                this 
                node    (1,:) = []
                policy  (1,1)   string  {mustBeMember(policy, ["Depth-First", "Width-First"])} = "Depth-First"
            end

            if this.isempty(), node = []; return; end

            if isempty(node), node = this.node_active; end

            switch policy
                case "Depth-First"
                    nodes = mStack();
                    for k = 1:numel(node.Children)
                        nodes.push(node.Children(k));
                    end
                    node = nodes.pop();
                case "Width-First"
                    nodes = mQueue();
                    for k = 1:numel(node.Children)
                        nodes.enqueue(node.Children(k));
                    end
                    node = nodes.dequeue();
                otherwise
                    node = [];
            end
        end

        % This function return if the active node is posterity of given
        % tree node or equal
        function tf = is_posterity_of(this, node)
            arguments
                this
                node (1,1)  matlab.ui.container.TreeNode
            end

            if node == this.node_active
                tf = true;
            else
                tf = false;

                node = node.Parent;
                while ~isempty(node) && isa(node, "matlab.ui.container.TreeNode")
                    if node == this.node_active, tf = true; break; end
                    node = node.Parent;
                end
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

            % update snapshot only
            this.update_snapshot_view(node);
        end

        % This function updates operation view
        % [Group] global control
        function update_operation_view(this)
            %% remove all style
            % can remove icon style only if Verison >= MATLAB R2024b
            if isMATLABReleaseOlderThan("R2024b")
                % remove previous font style on tree
                removeStyle(this.optree);

                fontStyle = uistyle("FontWeight", 'bold');
                addStyle(this.optree, fontStyle, "node", this.node_active);
                this.node_active.Icon = this.flagsrc;
                ndptr = this.get_node_before();
                while ~isempty(ndptr)
                    % add to previous node
                    addStyle(this.optree, fontStyle, "node", ndptr);
                    ndptr.Icon = this.readysrc;

                    % update
                    ndptr = this.get_node_before(ndptr);
                end
            else
                % remove all previous style on tree
                removeStyle(this.optree);

                % append bold style on current active branch
                fontStyle = uistyle("FontWeight", 'bold');
                iconStyle0 = uistyle("Icon",this.flagsrc, "IconAlignment",'leftmargin');
                iconStyle = uistyle("Icon",this.readysrc, "IconAlignment",'leftmargin');
                addStyle(this.optree, fontStyle, "node", this.node_active);
                addStyle(this.optree, iconStyle0, "node", this.node_active);
                ndptr = this.get_node_before();
                while ~isempty(ndptr)
                    % add to previous node
                    addStyle(this.optree, fontStyle, "node", ndptr);
                    addStyle(this.optree, iconStyle, "node", ndptr);

                    % update
                    ndptr = this.get_node_before(ndptr);
                end
            end

            %% expand older node
            pnode = this.get_node_before();
            if ~isempty(pnode), pnode.expand(); end

            %% update current selected node
            this.optree.SelectedNodes = this.node_active;
        end

        % This function updates source view
        % [Group] global control
        function update_source_view(this)
            % invoke storage agent scan the project

            % and repaint
        end

        % This function updates snapshot view
        % [Group] local control
        function update_snapshot_view(this, node)
            arguments
                this
                node    (1,:)   = []
            end

            % set default node as selected node
            if isempty(node), node = this.optree.SelectedNodes; end
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

    methods (Static)
        function r = reformat_time(opt_crop)
            frames = opt_crop.f;
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
end

