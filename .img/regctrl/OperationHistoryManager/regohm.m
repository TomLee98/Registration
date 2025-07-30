classdef regohm < handle
    %REGOHM This class defines "Operation History Manager" for reg3D app, 
    % which is the agent to handle image and operations resources

    properties(Constant, Hidden)
        OP_SKIP_OUTPUT = {constdef.OP_LOAD, constdef.OP_SEGMENT}
    end

    properties(Access = public, Dependent)
        ActiveNodeData  % ___/get, 1-by-1 struct with field {arg, cbot, data, optr}
        ActiveNodeTag   % ___/get, 1-by-1 string, could be empty
        CachePolicy     % ___/get, 1-by-1 string, could be "PERFORMANCE"/"RESOURCE"/"BALANCE"
        CacheSizeMEM    % set/get, 1-by-1 double, indicate the maximal memory cache size
        CacheSizeHDD    % set/get, 1-by-1 double, indicate the maximal hard drive cache size
        CurrentImagePtr % ___/get, 1-by-1 regmov, indicate the current activated image
    end

    properties(SetAccess = immutable, GetAccess = private)
        ufig            % 1-by-1 matlab.ui.Figure
        optree          % 1-by-1 matlab.ui.container.Tree
        optxt           % 1-by-1 matlab.ui.controller.TextArea
        opsmgr          % 1-by-1 storage_mananer object (component app)
        ctmenu          % 1-by-1 matlab.ui.container.ContextMenu
        flagsrc = string(fileparts(mfilename("fullpath"))).extractBefore(...
            "regctrl") + filesep + "sources" + filesep + "active.png"
        readysrc = string(fileparts(mfilename("fullpath"))).extractBefore(...
            "regctrl") + filesep + "sources" + filesep + "ready.png"
        calcsrc = string(fileparts(mfilename("fullpath"))).extractBefore(...
            "regctrl") + filesep + "sources" + filesep + "rtcalc.png"
    end

    properties(Access = private, Hidden)
        cache_size_mem  (1,1)   double  = 64                % memory cache size, unit as GB 
        cache_size_hdd  (1,1)   double  = 128               % hard drive cache size, unit as GB
        dptr            (1,1)   regmov = regmov.empty()     % current node related data pointer
        node_active     (1,1)                               % current active node (or tree)
        nodes_total_num (1,1)   double  = 0                 % number of total generated nodes in the tree
        optsty          (1,1)   string  = "PERFORMANCE"     % 1-by-1 string, indicate the manager running method
        storage_summary (1,1)   struct                      % struct with storage summary
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

        function r = get.CacheSizeMEM(this)
            r = this.cache_size_mem;
        end

        function set.CacheSizeMEM(this, r)
            arguments
                this
                r       (1,1)   double  {mustBeNonnegative}
            end

            aes = aesobj();
            eval(aes.decrypt(constdef.MEM_CACHE_KEY));
            if r <= MEM_CACHE_SIZE_MAX
                this.cache_size_mem = r;
            else
                throw(MException("regohm:invalidCacheSize", ...
                    "Please contact the administrator to allocate more cache."));
            end

            % update view
            this.update_storage_summary();
            this.update_storage_view([1,2]);
        end

        function r = get.CacheSizeHDD(this)
            r = this.cache_size_hdd;
        end

        function set.CacheSizeHDD(this, r)
            arguments
                this
                r       (1,1)   double  {mustBeNonnegative}
            end

            aes = aesobj();
            eval(aes.decrypt(constdef.HDD_BUFFER_KEY));
            if r <= HDD_BUFFER_SIZE_MAX
                this.cache_size_hdd = r;
            else
                throw(MException("regohm:invalidCacheSize", ...
                    "Please contact the administrator to allocate more cache."));
            end

            % update view
            this.update_storage_summary();
            this.update_storage_view([1,2]);
        end

        function r = get.CachePolicy(this)
            r = this.optsty;
        end

        function r = get.CurrentImagePtr(this)
            r = this.dptr;
        end
        
    end

    methods (Access = public)
        function this = regohm(ufig, op_tree, op_txt, op_smgr, ctxt_menu, cache_mem, cache_hdd, strategy)
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
                op_smgr     (1,1)   storage_manager
                ctxt_menu   (1,1)   matlab.ui.container.ContextMenu
                cache_mem   (1,1)   double
                cache_hdd   (1,1)   double
                strategy    (1,1)   string  {mustBeMember(strategy, ...
                                    ["PERFORMANCE", "RESOURCE", "BALANCE"])} = "PERFORMANCE"
            end

            this.ufig = ufig;
            this.optree = op_tree;
            this.optree.SelectionChangedFcn = @(~, event)this.sltchg_callback(event);
            this.optxt = op_txt;
            this.opsmgr = op_smgr;
            this.opsmgr.SetBoss(this);
            this.ctmenu = ctxt_menu;
            this.optsty = strategy;
            this.node_active = this.optree;

            % validate
            this.CacheSizeMEM = cache_mem;
            this.CacheSizeHDD = cache_hdd;
            % pass

            this.storage_summary.Overview = array2table([cache_mem, cache_hdd; ...
                                                         0,         0; ...
                                                         0,         0], ...
                "VariableNames",["MEM","HDD"], "RowNames",["Full", "Used", "UsedRatio"]);
            this.storage_summary.Details = table('Size',[0,3], ...
                'VariableTypes',{'string','double','string'}, ...
                'VariableNames',{'NodeTag', 'Used', 'Location'});
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

        function tf = isaddable(this, optr, args)
            % note that current operation generates data must be with less
            % or equal size than active node data
            % And there must be enough space to allow current image display
            arguments
                this 
                optr (1,1)  string  {mustBeMember(optr, ["CROP","LOAD","SEGMENT","REGISTER"])}
                args (1,:)  cell    = {}        % Name-Value pair, example: "Z",[5, 10]
            end

            % empty tree, must be addable
            if this.isempty() == true, tf = true; return; end

            % storage policy determines if new storage space needed
            is_storage = constdef.TRANSMISSION_STORAGE_TABLE{this.optsty, optr};
            if is_storage== false, tf = true; return; end

            % calculate storage free size
            mem_left = this.storage_summary.Overview.MEM(1) ...
                - this.storage_summary.Overview.MEM(2);     % GB
            hdd_left = this.storage_summary.Overview.HDD(1) ...
                - this.storage_summary.Overview.HDD(2);     % GB

            % find the nearest stored node and get the storage size
            node = this.find_nearest_prior_stored_node(this.node_active);
            dptr_tmp = node.NodeData.RSPoint.getData();
            alloc_st = (dptr_tmp.Bytes.mem + dptr_tmp.Bytes.map)/2^30;   % GB

            % use operator and arguments to estimate the storage using
            switch optr
                case constdef.OP_CROP
                    switch args{1}
                        case "XY"
                            alloc_new = prod(diff(args{2},1,1)+1) ...
                                /(dptr_tmp.MetaData.width*dptr_tmp.MetaData.height) ...
                                * alloc_st;
                        case "Z"
                            alloc_new = (diff(args{2},1)+1)/dptr_tmp.MetaData.slices ...
                                * alloc_st;
                        case "T"
                            alloc_new = numel(args{2})/dptr_tmp.MetaData.frames ...
                                * alloc_st;
                        otherwise
                    end

                    tf = alloc_new < max(mem_left, hdd_left);
                case constdef.OP_LOAD
                    % loading must be addable
                    tf = true;
                case constdef.OP_SEGMENT
                    % note that NuclearCtrlBot and some variables in the
                    % memory, and space using could be omitted
                    tf = true;
                case constdef.OP_REGISTER
                    % note that register need space as active node if
                    % policy is "performance"
                    tf = (alloc_st < max(mem_left, hdd_left));
                otherwise
            end
        end

        function AddNode(this, node, optr, args, varargin)
            % This function create new node at optree and store a regrspt
            % object in the node

            %% append node after active node
            node_new = this.create_new_node(node, optr, args, varargin);

            %% activate new node
            this.ActivateNode(node_new);

            %% update storage summary
            this.update_storage_summary();

            this.update_manage_view();
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

            %% update storage summary
            this.update_storage_summary();

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
               this.update_operation_view();
               this.update_snapshot_view();
            end
        end

        function MergeNodes(this, skflag)
            arguments
                this 
                skflag  (1,1)   logical = true
            end

            %% validate if data size and operation are compatible
            slt_nodes = this.optree.SelectedNodes;  % 
            for k = 1:numel(this.optree.SelectedNodes)-1
                nd_cur = this.optree.SelectedNodes(k).NodeData;
                nd_next = this.optree.SelectedNodes(k+1).NodeData;
                if any(nd_cur.RSPoint.ImageDim ~= nd_next.RSPoint.ImageDim)
                    % ERROR: display at anytime
                    uialert(this.ufig, "Nodes must be with equal data size.", ...
                        "Error", "Icon","error");
                    return;
                end

                opt_cur = nd_cur.RSPoint.Arguments.(constdef.OP_REGISTER);
                opt_next = nd_next.RSPoint.Arguments.(constdef.OP_REGISTER);

                if (opt_cur.Algorithm==opt_next.Algorithm) ...
                        || (ismember(opt_cur.Algorithm, ["TCREG","OCREG","MANREG"]) ...
                        && ismember(opt_next.Algorithm, ["TCREG","OCREG","MANREG"]))
                    % continue
                else
                    % ERROR: display at anytime
                    uialert(this.ufig, "Nodes must be with compatible registration type.", ...
                        "Error", "Icon","error");
                    return;
                end
            end

            %% create combined node
            % generate nodes line
            [~, nodes] = find_nearest_prior_stored_node(this, ...
                this.optree.SelectedNodes(1));
            % real time calculation chain
            node = nodes.pop();
            rpt = node.NodeData.RSPoint;
            dptr_tmp = rpt.getData();     % get storage, will be thrown after merging
            while ~isempty(nodes)
                node = nodes.pop();
                rpt = node.NodeData.RSPoint;

                % update current dptr, source only one copy
                dptr_tmp = rpt.getData(dptr_tmp);
            end

            % generate new data storage
            img_comb = dptr_tmp.copy();     % current data deep copy

            % generated variables
            args_comb = this.optree.SelectedNodes(1).NodeData.RSPoint.Arguments;
            tfs_comb = this.optree.SelectedNodes(1).NodeData.RSPoint.Transformation;
            others_comb = this.optree.SelectedNodes(1).NodeData.RSPoint.Others;
            for k = 2:numel(this.optree.SelectedNodes)
                % combined transformation
                vloc = cellfun(@(x)isempty(x), tfs_comb, "UniformOutput",true);
                tfs_k = this.optree.SelectedNodes(k).NodeData.RSPoint.Transformation;
                ploc = cellfun(@(x)~isempty(x), tfs_k, "UniformOutput",true);
                tfs_comb(vloc&ploc) = tfs_k(vloc&ploc);

                % modify others
                others_k = this.optree.SelectedNodes(k).NodeData.RSPoint.Others;
                others_comb.fixdef.Sampling(2) = ...
                    string(unique([others_comb.fixdef.Sampling(2), others_k.fixdef.Sampling(2)])).join(",");
            end
            ploc = cellfun(@(x)~isempty(x), tfs_comb, "UniformOutput",true);
            optf.f = find(ploc);
            r = regohm.reformat_time(optf);
            others_comb.frames_reg(1) = string(r).join(":").join(",");

            % add new node and activate
            this.AddNode(this.optree, constdef.OP_REGISTER, args_comb, "Data",img_comb, ...
                "Transform",tfs_comb, "Others", others_comb);

            if ~skflag
                % remove nodes
                for k = 1:numel(slt_nodes)
                    node = slt_nodes(k);
                    this.DelNode(node);
                end
            end

            %% update appearance
            this.update_storage_summary();
            this.update_manage_view();

        end

        function ActivatePreviousNode(this)
            node = get_node_before(this);

            if ~isempty(node)
                this.ActivateNode(node); 
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
            node = this.optree.SelectedNodes(end);
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

        % This function find node with tag <node> and modify the storage
        function TrySwitchStorage(this, tag, storage)
            arguments
                this 
                tag     (1,1)   string
                storage (1,1)   string  {mustBeMember(storage, ["MEM","HDD"])}
            end

            % use depth-first traversal to spread data onto disk
            st = mStack();
            st.push(this.optree);   % push the tree root
            while ~isempty(st)
                % pop the node
                node = st.pop();

                if isa(node, "matlab.ui.container.TreeNode")
                    if node.Tag == tag
                        udlg = uiprogressdlg(this.ufig, "Indeterminate","on", ...
                            "Title","Message", "Icon","info");

                        switch storage
                            case "MEM"
                                udlg.Message = "Gathering (HDD/SSD->Memory) ...";
                                regohm.node_storage_switch(node, "gather");
                            case "HDD"
                                udlg.Message = "Spreading (Memory->HDD/SSD) ...";
                                regohm.node_storage_switch(node, "spread");
                            otherwise
                        end

                        delete(udlg);

                        % tag is unique, just break
                        break;
                    end
                end

                nodes = node.Children;     % nodes
                for k = 1:numel(nodes), st.push(nodes(k)); end
            end

            this.update_storage_summary();

            % update view
            this.update_storage_view([1,2]);
            this.update_snapshot_view();
        end

        function Save(this, file)
            arguments
                this
                file    (1,1)   string
            end

            % This function uses struct to save the operation tree
            % construct tree-struct dynamically
            [~, ~, ext] = fileparts(file);
            switch ext
                case constdef.PROJECT_FILE_EXT
                    DS_ = struct("cache_size_mem",  this.cache_size_mem, ...
                                 "cache_size_hdd",  this.cache_size_hdd, ...
                                 "dptr",            this.dptr, ...
                                 "optsty",          this.optsty, ...
                                 "tag_active",      this.node_active.Tag, ...
                                 "nodes_total_num", this.nodes_total_num, ...
                                 "storage_summary", this.storage_summary);

                    op_tree = opTree(this.optree, DS_);
                    save(file, "op_tree", '-mat','-append', '-nocompression');
                otherwise
                    throw(MException("regohm:invalidFileFormat", ...
                        "File extension must be %s.", constdef.PROJECT_FILE_EXT));
            end
        end

        function Load(this, file)
            arguments
                this
                file    (1,1)   string  {mustBeFile}
            end

            % This function reconstruct regohm from given file
            [~, ~, ext] = fileparts(file);
            switch ext
                case constdef.PROJECT_FILE_EXT
                    load(file, '-mat', "op_tree");
                    op_tree.RestoreUITree(this.optree, this.ctmenu);
                    DS_ = op_tree.OhmgrProperties;
                    this.cache_size_mem = DS_.cache_size_mem;
                    this.cache_size_hdd = DS_.cache_size_hdd;
                    this.dptr = DS_.dptr;
                    this.optsty = DS_.optsty;
                    this.nodes_total_num = DS_.nodes_total_num;
                    this.storage_summary = DS_.storage_summary;
                    node = this.find(DS_.tag_active);
                    if ~isempty(node)
                        this.ActivateNode(node);
                    else
                        this.node_active = this.optree;
                    end
                otherwise
                    throw(MException("regohm:invalidFileFormat", ...
                        "File extension must be %s.", constdef.PROJECT_FILE_EXT));
            end
        end
    end

    methods (Access = private)

        % This function creates a new node as children of active node
        % note that this function also determines storage
        function node = create_new_node(this, node, optr, args, vars)

            %% create restore point object by cache policy
            is_storage = constdef.TRANSMISSION_STORAGE_TABLE{this.optsty, optr};

            if is_storage == false
                % drop data node and replace by empty
                for n = 1:2:numel(vars)
                    if isequal("Data", vars{n})
                        % drop, handle is only hold by Register now
                        % replace with a place holder
                        tmprv = vars{n+1};
                        vars{n+1} = regmov.make_placeholder_as(tmprv);
                        break;
                    end
                end
            end

            % construct restore point
            % note that the data must be on disk
            rs_node = regrspt(optr, args, vars{:});

            node_data = struct("Operation",  [], ...
                               "Properties", [], ...
                               "RSPoint",    rs_node, ...
                               "Time",       datetime("now", "Format","MM/dd HH:mm:ss"));

            switch optr
                case constdef.OP_LOAD
                    node_text = "load";
                    % must be on hard-drive
                    node_data.Properties.Location = tmprv.Location;
                    node_data.Properties.Dims = [tmprv.MetaData.width, ...
                        tmprv.MetaData.height, tmprv.MetaData.channels, ...
                        tmprv.MetaData.slices, tmprv.MetaData.frames];
                    node_data.Properties.Size = round(tmprv.Bytes.map/2^30, 1);
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
                    node_text = sprintf("register(%s)", st.Mode);
                    switch st.Algorithm
                        case {"TCREG", "OCREG"}
                            switch st.Mode
                                case "global"
                                    node_data.Properties = struct("Algorithm",      st.Algorithm, ...
                                                                  "RegisterFrames", rs_node.Others.frames_reg(1), ...
                                                                  "Transformation", st.Options.TformType, ...
                                                                  "CoarseRegister", st.Options.CoarseAlg, ...
                                                                  "MaxStep",        string(st.Options.MaxStep), ...
                                                                  "MinStep",        string(st.Options.MinStep));
                                case "local"
                                    node_data.Properties = struct("Algorithm",          st.SubAlgorithm, ...
                                                                  "RegisterFrames", rs_node.Others.frames_reg(1), ...
                                                                  "Compensation",       string(st.Options.ImageRehist), ...
                                                                  "CompensationGrade",  string(st.Options.RepAcc));
                                otherwise
                            end
                        case "MANREG"
                            node_data.Properties = struct("Transformation", st.Options.TformType, ...
                                                          "RegisterFrames", rs_node.Others.frames_reg(2), ...
                                                          "Degree",         string(st.Options.Degree), ...
                                                          "ProjectedView",  st.Options.DView, ...
                                                          "Isometric",      string(st.Options.Isometric));
                        case "LTREG"
                            node_data.Properties = struct("Keyframes",      string(reshape(st.Options.Keyframes,1,[])).join(","), ...
                                                          "RegisterFrames", rs_node.Others.frames_reg(1), ...
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

            %% add node to active node or given node
            if isempty(node)
                node = uitreenode(this.node_active, "Text",node_text, ...
                    "NodeData",node_data, "ContextMenu",this.ctmenu);
            else
                if isvalid(node) && (isa(node, "matlab.ui.container.TreeNode") ...
                        || isa(node, "matlab.ui.container.Tree"))
                    node = uitreenode(node, "Text",node_text, ...
                        "NodeData",node_data, "ContextMenu",this.ctmenu);
                end
            end
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
                % generate nodes line
                [~, nodes] = find_nearest_prior_stored_node(this, node);

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

                    % update current dptr, source only one copy
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

        % This function finds the node which has Tag equals to given tag
        function node = find(this, tag)
            arguments
                this
                tag     (1,1)   string
            end

            node = matlab.graphics.GraphicsPlaceholder.empty();
            
            nodes = mQueue();
            nodes.enqueue(this.optree);
            while ~isempty(nodes)
                nd = nodes.dequeue();
                if isa(nd, "matlab.ui.container.TreeNode")
                    if nd.Tag == tag, node = nd; return; end
                end
                
                for k = 1:numel(nd.Children)
                    nodes.enqueue(nd.Children(k));
                end
            end
        end

        % This function auto moves the active node if it is on the branch
        % which will be deleted
        function automove_active_node(this, node)
            if this.is_posterity_of(node)   % active_node is posterity of node?
                % priority: brother > parent
                if numel(node.Parent.Children) > 1
                    brothers = setdiff(node.Parent.Children, node);
                    this.ActivateNode(brothers(1));
                else
                    if isa(node.Parent, "matlab.ui.container.TreeNode")
                        % activate the nearest
                        this.ActivateNode(node.Parent);
                    else
                        % only one node, deactivate
                        this.deactivate();
                    end
                end
            end
        end

        % This function update manager appearance accords to current 
        % tree status
        function update_manage_view(this)
            % global control
            this.update_operation_view();
            this.update_storage_view();

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
                if isvalid(this.node_active) && ...
                        isa(this.node_active.Parent, "matlab.ui.container.TreeNode")
                    node = this.node_active.Parent;
                else
                    node = [];
                end
            else
                if isvalid(node) && isa(node.Parent, "matlab.ui.container.TreeNode")
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

                node_tmp = this.node_active.Parent;
                while ~isempty(node_tmp) && isa(node_tmp, "matlab.ui.container.TreeNode")
                    if node == node_tmp, tf = true; break; end
                    node_tmp = node_tmp.Parent;
                end
            end
        end

        % find the shortest calculation path, equals to find the
        % nearest storage which is of a higher seniority in the
        % family hierarchy
        function [node, nodes] = find_nearest_prior_stored_node(this, node)
            arguments
                this 
                node (1,1)  matlab.ui.container.TreeNode
            end

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

            node = nodes.top();
        end

        % This is select changed callback binding on tree dynamically
        function sltchg_callback(this, event)
            if ~isempty(event.SelectedNodes)
                % get node data (select the last node)
                node = event.SelectedNodes(end);

                % update snapshot only
                this.update_snapshot_view(node);
            end
        end

        % This function updates operation view
        % [Group] global control
        function update_operation_view(this)

            if this.isempty(), return; end

            %% remove all style
            % can remove icon style only if Version >= MATLAB R2024b
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
            if isvalid(this.node_active)
                this.optree.SelectedNodes = this.node_active;
            else
                this.optree.SelectedNodes = [];
            end
        end

        % This function updates source view
        % [Group] global control
        function update_storage_view(this, comp)
            arguments
                this 
                comp (1,:)  double  {mustBeMember(comp, [1,2,3])} = [1,2,3]
            end
            % invoke storage manager agent repaint
            this.opsmgr.UpdateView(this.storage_summary, comp);
        end

        % This function updates snapshot view
        % [Group] local control
        function update_snapshot_view(this, node)
            arguments
                this
                node    (1,:)   = []
            end

            % set default node as selected node
            if isempty(node)
                if ~isempty(this.optree.SelectedNodes)
                    node = this.optree.SelectedNodes(end);
                else
                    if ~this.isempty()
                        node = this.node_active;
                    else
                        node = [];  % no node left
                    end
                end
            end

            if isempty(node)
                txt = "";
            else
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
            end

            %% update text field to display the snapshot
            this.optxt.Value = txt;
        end

        % This function visit all tree and summary the storage view
        function update_storage_summary(this)
            if ~this.isempty()
                % use depth-first traversal to spread data onto disk
                % TODO: if there is too many nodes, use incremental update
                NodeTag = string([]);
                Used = [];
                Location = string([]);
                MEM_Used = 0;
                HDD_Used = 0;

                st = mStack();
                st.push(this.optree);   % push the tree root
                while ~isempty(st)
                    % pop the node
                    node = st.pop();

                    if isa(node, "matlab.ui.container.TreeNode")
                        nd = node.NodeData;
                        % get node data summary
                        NodeTag = [NodeTag; node.Tag]; %#ok<AGROW>
                        dptr_tmp = nd.RSPoint.getData();
                        switch dptr_tmp.Location
                            case "Memory"
                                Location = [Location; "MEM"]; %#ok<AGROW>
                                Used = [Used; dptr_tmp.Bytes.mem/2^30]; %#ok<AGROW>
                                MEM_Used = MEM_Used + dptr_tmp.Bytes.mem;
                            otherwise
                                Location = [Location; "HDD"]; %#ok<AGROW>
                                Used = [Used; dptr_tmp.Bytes.map/2^30]; %#ok<AGROW>
                                HDD_Used = HDD_Used + dptr_tmp.Bytes.map;
                        end
                    end

                    nodes = node.Children;     % nodes
                    for k = 1:numel(nodes), st.push(nodes(k)); end
                end

                % summary to table
                this.storage_summary.Overview = array2table([this.cache_size_mem, this.cache_size_hdd; ...
                                                             MEM_Used/2^30,       HDD_Used/2^30; ...
                                                             MEM_Used/(this.cache_size_mem*2^30), HDD_Used/(this.cache_size_hdd*2^30)], ...
                    "VariableNames",["MEM","HDD"], "RowNames",["Full", "Used", "UsedRatio"]);
                this.storage_summary.Details = table(NodeTag, Used, Location);
            else
                % generate empty table
                this.storage_summary.Overview = array2table([this.cache_size_mem, this.cache_size_hdd; ...
                                                             0,                   0; ...
                                                             0,                   0], ...
                    "VariableNames",["MEM","HDD"], "RowNames",["Full", "Used", "UsedRatio"]);
                this.storage_summary.Details = table('Size',[0,3], ...
                    'VariableTypes',{'string','double','string'}, ...
                    'VariableNames',{'NodeTag', 'Used', 'Location'});
            end
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

        function node_storage_switch(node, optr)
            arguments
                node    (1,1)   matlab.ui.container.TreeNode
                optr    (1,1)   string  {mustBeMember(optr, ["spread", "gather"])}
            end

            nd = node.NodeData;

            switch optr
                case "spread"
                    nd.RSPoint.IsOnHardDrive = true;    % spread data
                case "gather"
                    nd.RSPoint.IsOnHardDrive = false;    % spread data
                otherwise
            end
        end
    end
end

