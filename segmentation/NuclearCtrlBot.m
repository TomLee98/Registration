classdef NuclearCtrlBot < handle
    %VIEWER This class is nuclear viewer object definition
    
    properties(Constant, Hidden)
        STATUS_SUCCESS      =   0
        STATUS_NO_HANDLE    =   -1
        REMOVED_OBJ         =   -1
        STATUS_DUPLICATE_OBJ=   -2
        
        VALID_FONT_COLOR    =   [1, 1, 1]       % valid object font color
        NORMAL_CIRCLE_COLOR =   [.5, .5, .5]
    end

    properties(SetAccess=private, GetAccess=private, Hidden)
        caller                              % app caller, must be ReTiNA object
        ax_caller                           % axes from caller, must be uiaxes
        nuclears                            % NuclearGroup object, for objects record
        hobj        (:,1)   cell            % handle of circle objects
        hrange      (:,1)           = []    % handle of display range (rectangle)
        cmap                                % string, indicate the color map
        disp_flag   (1,1)   string  = "on"  % string, indicate ROIs display state, "on"/"off"
        disp_range  (:,4)           = []    % 0/1-by-4 positive integer double, as [x0, y0, w, h]
    end

    properties(SetAccess=private, GetAccess=private)
        nuid            % nuclear indicator definition(identity)
        hlid            % 1-by-s double, the highlight id array
        volopts         % the image options of vol
        center          % nuclear indicator center
        radius          % nuclear indicator radius
        color           % the color map for viewer, m-by-3 matrix
        rois            % ROIsInfo struct
    end

    properties(Access=public, Dependent)
        DispFlag
        DispRange
        HighlightID
        IsMaskExist
        MaxN
        Parent
        ROIsInfo
    end

    methods
        function r = get.DispFlag(this)
            r = this.disp_flag;
        end

        function set.DispFlag(this, df)
            arguments
                this 
                df      (1,1)   string  {mustBeMember(df, ["on","off"])}
            end

            this.disp_flag = df;
        end

        function r = get.DispRange(this)
            r = this.disp_range;
        end

        function set.DispRange(this, r)
            arguments
                this
                r       (:,4)   double  {mustBeNonnegative, mustBeInteger} = []
            end

            this.disp_range = r;
        end

        function r = get.HighlightID(this)
            r = this.hlid;
        end

        function set.HighlightID(this, r)
            arguments
                this 
                r   (1,:)   double  {mustBePositive, mustBeInteger}
            end
            
            if ~isempty(r)
                if max(r) <= this.nuclears.obj_max_exist
                    this.hlid = r;
                else
                    throw(MException("NuclearCtrlBot:invalidID", ...
                        "The given id out of boundary."));
                end
            else
                this.hlid = r;
            end
        end

        function r = get.IsMaskExist(this)
            r = ~isempty(this.nuclears);
        end

        function r = get.Parent(this)
            r = this.caller;
        end

        function r = get.ROIsInfo(this)
            r = this.rois;
        end

        function r= get.MaxN(this)
            if ~isempty(this.nuclears.obj_max_exist)
                r = this.nuclears.obj_max_exist;
            else
                r = zeros(1, size(this.nuclears.obj_max_exist, 2));
            end
        end
    end
    
    methods(Access = public)
        function this = NuclearCtrlBot(res, vopt, cmap)
            %NUVIEWER A Constructor
            arguments
                res     (1,1)   struct
                vopt    (1,6)   table
                cmap    (1,1)   string {mustBeMember(cmap, ...
                    ["default","hsv","jet","parula","turbo"])} = "default"
            end

            this.center = res.Center;
            this.radius = res.Radius;
            this.nuid = res.Nid;
            this.hlid = zeros(1, 0);
            this.rois = GenROI(this.nuid, this.center, this.radius);
            this.volopts = vopt;
            this.cmap = cmap;
            this.nuclears = NuclearGroup(this.nuid);
            this.gen_colormap();
        end

        function status = Display(this, zidx)
            arguments
                this
                zidx (1,1) double {mustBeNonnegative, mustBeInteger} = 1
            end

            %DISPLAY Disp the objects on viewer ax
            if isempty(this.ax_caller) || ~isvalid(this.ax_caller)
                status = this.STATUS_NO_HANDLE;
                return;
            end

            % skip
            if this.disp_flag ~= "on"
                status = this.STATUS_SUCCESS;
                return;
            end

            % create display range if needed
            if ~isempty(this.disp_range)
                if isempty(this.hrange) || ~isvalid(this.hrange)
                    this.hrange = images.roi.Rectangle(this.ax_caller, ...
                        "Position",this.disp_range, "Color","y", ...
                        "FaceAlpha",0, "StripeColor",[0.1, 0.1, 0.1], "Deletable",false, ...
                        "InteractionsAllowed","none", "LineWidth",1, "MarkerSize",0.1);
                else
                    % keep exist rectangle ROI
                end
            else
                if ~isempty(this.hrange) && isvalid(this.hrange)
                    delete(this.hrange);
                end

                this.hrange = [];
            end

            id_z = this.nuid{zidx};
            if ~isempty(id_z)
                h_obj = cell(numel(id_z),1);
                % draw circle at the plane
                for cir_index = 1:numel(id_z)
                    cir_id = id_z(cir_index);

                    % skip the removed object
                    if cir_id == this.REMOVED_OBJ, continue; end

                    if ~isnan(cir_id)
                        if ismember(cir_id, this.hlid)
                            cc = this.color{zidx}(cir_index,:);
                        else
                            cc = this.NORMAL_CIRCLE_COLOR;   % deep gray
                        end
                        % validate if circle in range (consider center for simple)
                        if isempty(this.disp_range)
                            create_new_object = true;
                        else
                            x0 = this.center{zidx}(cir_index, 1);
                            y0 = this.center{zidx}(cir_index, 2);
                            create_new_object = (x0 > this.disp_range(1) ...
                                && x0 < this.disp_range(1)+this.disp_range(3)) && ...
                                (y0 > this.disp_range(2) ...
                                && y0 < this.disp_range(2)+this.disp_range(4));
                        end

                        if create_new_object == true
                            self = images.roi.Circle(this.ax_caller,...
                                "Center",this.center{zidx}(cir_index,:),...
                                "Radius",this.radius{zidx}(cir_index), ...
                                "Color",cc, ...
                                "LabelTextColor", this.VALID_FONT_COLOR, ...
                                "LineWidth",1.5, "MarkerSize",1, "Label",string(cir_id),...
                                "LabelVisible","on", "FaceAlpha",0, "LabelAlpha",0,...
                                "Deletable",true, "UserData",[cir_index, zidx]);

                            % binding listener(will be auto removed if the roi was deleted)
                            addlistener(self, ...
                                'ROIMoved', @(src,~)this.roi_moved(src,[],self));

                            % binding menu: 修改ID
                            uimenu(self.ContextMenu, ...
                                "Text","修改ID",...
                                "MenuSelectedFcn",@(~,~)this.roi_modify_id([],[],self), ...
                                "Tag","IDROIContextMenuModify");

                            % changed the default context menu: delete
                            self.ContextMenu.Children(end).MenuSelectedFcn ...
                                = @(~,~)this.roi_delete([],[],self);

                            h_obj{cir_index} = self;
                        end
                    end
                end
            else
                h_obj = {};
            end

            this.hobj = h_obj;
            status = this.STATUS_SUCCESS;
        end

        function LoadMaskFrom(this, file)
            % This function loads mask from file
            arguments
                this 
                file    (1,1)   string  {mustBeFile}
            end

            try
                % load mat file
                load(file, "-mat", "Center", "Radius", "Nid", "VolOpts", "CMap");
            catch ME
                throw(ME);
            end

            % update
            this.center = Center;
            this.radius = Radius;
            this.nuid = Nid;
            this.hlid = zeros(1, 0);
            this.rois = GenROI(this.nuid, this.center, this.radius);
            
            this.volopts = VolOpts;
            this.cmap = CMap;

            this.disp_flag = "on";
            this.nuclears = NuclearGroup(this.nuid);
            this.gen_colormap();
        end

        function SaveMaskTo(this, file)
            % This function saves mask to file
            arguments
                this 
                file (1,1)  string
            end

            Center = this.center;
            Radius = this.radius;
            Nid = this.nuid;
            VolOpts = this.volopts;
            CMap = this.cmap;

            try
                % save as mat file
                save(file, "Center","Radius","Nid","VolOpts","CMap","-mat");
            catch ME
                throw(ME);
            end
        end

        function mask = GenMaskVol(this)
            % This function uses rois to create volume mask
            vsize = [this.volopts.height, this.volopts.width, this.volopts.slices];
            vroi = ROIP2V(this.rois, vsize);

            % transform as volume coding
            Vm = zeros(vsize, "uint16");
            for k = 1:numel(vroi), Vm(vroi{k}) = k; end

            % transform to sparse coding
            row_mask = [];
            col_mask = [];
            % transform the format: [comp1, comp2, ..., comp_n]
            % where comp_k is column array with 1 indicate voxels
            d = unique(Vm);
            d(d==0) = [];
            d(isnan(d)) = [];

            for k = 1:numel(d)
                vd_loc = find(Vm == d(k));   % F order, linear sub index
                row_mask = [row_mask; vd_loc]; %#ok<AGROW>
                col_mask = [col_mask; k*ones(size(vd_loc))]; %#ok<AGROW>
            end
            vol_mask = ones(size(row_mask));

            % generate sparse matrix
            mask = sparse(row_mask, col_mask, vol_mask, numel(Vm), numel(d));
        end

        function status = Clear(this)
            if ~isempty(this.hrange) && isvalid(this.hrange)
                delete(this.hrange);
            end
            this.hrange = [];

            if ~isempty(this.hobj)
                for k = 1:numel(this.hobj)
                    delete(this.hobj{k}); 
                end
                status = this.STATUS_SUCCESS;
            else
                status = this.STATUS_NO_HANDLE;
            end
        end

        function status = AddObject(this, zidx, cir_id, cir_info)
            % this function add an object to data
            arguments
                this
                zidx (1,1) double {mustBePositive, mustBeInteger}
                cir_id (1,1) double % could be nan for uncertainty
                cir_info (1,3) double {mustBePositive}  % [cx,cy,r]
            end

            if zidx > this.volopts.slices
                throw(MException("NuViewer:invalidSliceIndex", ...
                    "Slice index is out of range."));
            end

            if ~isnan(cir_id) && ismember(cir_id, this.nuid{zidx})
                status = this.STATUS_DUPLICATE_OBJ;
                return;
            else
                % add new record
                this.nuclears.AddObj(cir_id);
            end

            cir_index = numel(this.nuid{zidx});
            this.nuid{zidx}(cir_index+1,1) = cir_id;
            this.center{zidx}(cir_index+1, :) = cir_info(1:2);
            this.radius{zidx}(cir_index+1,1) = cir_info(3);

            this.rois = GenROI(this.nuid, this.center, this.radius);

            this.gen_colormap();

            status = this.STATUS_SUCCESS;
        end

        function status = Refresh(this)
            % this function reorder objects identities
            % access each object and modify the identities
            
            % apply the refreshed identity
            [this.nuid, this.center, this.radius] ...
                = this.nuclears.ApplyRefreshTo(this.nuid, ...
                this.center, this.radius, true);

            % refresh ROIs
            this.rois = GenROI(this.nuid, this.center, this.radius);

            % regenerate the color map
            this.gen_colormap();

            status = this.STATUS_SUCCESS;
        end

        function delete(this)
            Clear(this);

            % delete(this);
        end
    end

    methods(Access = ?CellSegmentor, Hidden)
        function SetHandle(this, caller)
            arguments
                this
                caller  (1,1)   Register
            end

            this.caller = caller;
            this.ax_caller = caller.UIAxes_Slice;
        end
    end

    methods(Access = private)
        function gen_colormap(this)
            % generate color map, colors number as same as valid object
            % number, store as center

            % get the total objects number
            total_n = this.nuclears.obj_num_tot;
            
            % generate the color map
            if this.cmap == "default"
                color_ = parula(total_n);
            else
                cfun = str2func(this.cmap);
                color_ = cfun(total_n);
            end

            % map the colors to objects
            this.color = cell(size(this.nuid));
            for zidx = 1:numel(this.nuid)
                id_z = this.nuid{zidx};
                for cidx = 1:numel(id_z)
                    id_z_c = id_z(cidx);
                    if ~isnan(id_z_c) && (id_z_c~=this.REMOVED_OBJ)
                        [~, id_loc] = this.nuclears.IsMember(id_z_c);
                        this.color{zidx}(cidx, :) ...
                            = color_(id_loc, :);
                    end
                end
            end
        end

        function status = move_object(this, zidx, cir_index, cir_info)
            % this function modify an object in data already
            arguments
                this
                zidx        (1,1) double {mustBePositive, mustBeInteger}
                cir_index   (1,1) double {mustBePositive, mustBeInteger}
                cir_info    (1,3) double {mustBePositive}  % [cx,cy,r]
            end
            if zidx > this.volopts.slices
                throw(MException("NuViewer:invalidSliceIndex", ...
                    "Slice index is out of range."));
            end
            if cir_index > numel(this.nuid{zidx})
                throw(MException("NuViewer:invalidObjectIndex", ...
                    "Object index is out of range."));
            end

            % find the object
            % identity: this.id{zidx}(cir_index)
            this.center{zidx}(cir_index, :) = cir_info(1:2);
            this.radius{zidx}(cir_index, :) = cir_info(3);

            status = this.STATUS_SUCCESS;
        end

        function status = recall_object(this, zidx, cir_index, cir_id_old, cir_id_new)
            % this function rename an object in data already
            arguments
                this
                zidx        (1,1) double {mustBePositive, mustBeInteger}
                cir_index   (1,1) double  {mustBeNonNan} % empty for not needed: cir_id_old is not nan
                cir_id_old  (1,1) double
                cir_id_new  (1,1) double
            end
            if zidx > this.volopts.slices
                throw(MException("NuViewer:invalidSliceIndex", ...
                    "Slice index is out of range."));
            end
            if isnan(cir_id_old) && isnan(cir_id_new)
                status = this.STATUS_SUCCESS;
                return;
            end
            if isempty(cir_index)
                if isnan(cir_id_old)
                    nanbi = isnan(this.nuid{zidx});
                    if sum(nanbi) > 1
                        % can't make sure the changed object
                        throw(MException("NuViewer:tooManyNaNIndex", ...
                            "Can't change label without index at this plane."));
                    else
                        % set nan -> val
                        this.nuid{zidx}(nanbi) = cir_id_new;
                        status = this.nuclears.NaN2Val(cir_id_new);
                        return;
                    end
                else
                    cir_index = (this.nuid{zidx}==cir_id_old);
                    if sum(cir_index) ~= 1
                        throw(MException("NuViewer:idNotFound", ...
                            "Identity not found."));
                    else
                        this.nuid{zidx}(cir_index) = cir_id_new;
                        status = this.nuclears.Val2Val(cir_id_old, cir_id_new);
                        return;
                    end
                end
            else
                if cir_index > numel(this.nuid{zidx})
                    throw(MException("NuViewer:invalidObjectIndex", ...
                        "Object index is out of range."));
                end
                cir_id_old = this.nuid{zidx}(cir_index);
                this.nuid{zidx}(cir_index,1) = cir_id_new;
                status = this.nuclears.Val2Val(cir_id_old, cir_id_new);
                return;
            end
        end

        function status = remove_object(this, zidx, cir_index)
            % this function remove an object from data
            arguments
                this
                zidx        (1,1) double {mustBePositive, mustBeInteger}
                cir_index   (1,1) double {mustBePositive, mustBeInteger}
            end

            if zidx > this.volopts.slices
                throw(MException("NuViewer:invalidSliceIndex", ...
                    "Slice index is out of range."));
            end
            if cir_index > numel(this.nuid{zidx})
                throw(MException("NuViewer:invalidObjectIndex", ...
                    "Object index is out of range."));
            end

            % marked the id set
            cir_id_old = this.nuid{zidx}(cir_index);
            status = this.nuclears.DelObj(cir_id_old);

            if status == this.STATUS_SUCCESS
                this.nuid{zidx}(cir_index,1) = this.REMOVED_OBJ;
            end

        end

        function roi_moved(this, src, ~, self)
            % call Move Object
            cidx = self.UserData(1);
            zidx = self.UserData(2);    % [circle_index, z_index]
            cinfo = [src.Center, src.Radius];

            this.move_object(zidx, cidx, cinfo);
        end

        function roi_modify_id(this, ~, ~, self)
            % generate a message box for new label acquirement
            cidx = self.UserData(1);
            zidx = self.UserData(2);    % [circle_index, z_index]
            id_cur = this.nuid{zidx}(cidx);

            prompt = {'Enter new identity:'};
            dlgtitle = 'rename';
            fieldsize = [1 45];
            definput = {num2str(id_cur)};

            answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
            if isempty(answer), return; end

            id_new = round(str2double(answer{1}));

            if ~isnan(id_new) && id_new >= 1
                this.recall_object(zidx, cidx, id_cur, id_new);
                % update this color
                this.gen_colormap();
                self.Label =  string(id_new);
                self.Color = this.color{zidx}(cidx, :);
            else
                warning("NuclearCtrlBot:invalidIdentity", ...
                    "Identity must be a positive integer.");
            end
        end

        function roi_delete(this, ~, ~, self)
            % remove the record stored in this
            cidx = self.UserData(1);
            zidx = self.UserData(2);    % [circle_index, z_index]
            this.remove_object(zidx, cidx);

            % remove the ROI object
            delete(self);
        end

        function clean_stack(this, roi_z)
            % Z ROI format, crop
            this.nuid = this.nuid(roi_z);
            this.center = this.center(roi_z);
            this.radius = this.radius(roi_z);
            this.color = this.color(roi_z);
        end
    end
end
