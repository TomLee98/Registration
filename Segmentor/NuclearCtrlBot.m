classdef NuclearCtrlBot < handle
    %VIEWER This class is nuclear viewer object definition
    
    properties(Constant, Hidden)
        STATUS_SUCCESS      =   0
        STATUS_NO_HANDLE    =   -1
        REMOVED_OBJ         =   -1
        STATUS_DUPLICATE_OBJ=   -2
        
        VALID_FONT_COLOR         =   [1, 1, 1]   % valid object font color
    end

    properties(SetAccess=private, GetAccess=private, Hidden)
        caller          % app caller, must be ReTiNA object
        ax_caller       % axes from caller, must be uiaxes
        nuclears        % NuclearGroup object, for objects record
        hobj            % handle of circle objects
        cmap            % string, indicate the color map
    end

    properties(SetAccess=private, GetAccess=private)
        nuid            % nuclear indicator definition(identity)
        nlbl            % nuclear labels
        volopts         % the image options of vol
        center          % nuclear indicator center
        radius          % nuclear indicator radius
        color           % the color map for viewer, m-by-3 matrix
        rois            % ROIsInfo struct
    end

    properties(Access=public, Dependent)
        MaxN
        ROIsInfo
        ROIsLabel
        Parent
    end

    methods
        function r = get.Parent(this)
            r = this.caller;
        end

        function r = get.ROIsInfo(this)
            r = this.rois;
        end

        function r= get.MaxN(this)
            r = this.nuclears.obj_max_exist;
        end

        function r = get.ROIsLabel(this)
            ids = reshape(unique(cell2mat(this.nuid), "rows", "sorted"),1,[]);
            ids(isnan(ids)) = [];
            ids(ids==this.REMOVED_OBJ) = [];

            r = strings(size(ids));
            for k = 1:numel(ids)
                r(k) = this.getLabel(ids(k));
            end
        end
    end
    
    methods(Access = public)
        function this = NuclearCtrlBot(nsper, vopt, cmap)
            %NUVIEWER A Constructor
            arguments
                nsper   (1,1)   NuclearSplitter
                vopt    (1,10)  table
                cmap    (1,1)   string {mustBeMember(cmap, ...
                    ["default","hsv","jet","parula","turbo"])} = "default"
            end

            this.center = nsper.Center;
            this.radius = nsper.Radius;
            this.nuid = nsper.Nid;
            this.nlbl = cellfun(@(x)string(x),this.nuid,"UniformOutput",false);
            this.rois = GenROI(this.nuid, this.center, this.radius);
            this.caller = nsper.Parent;
            this.ax_caller = this.caller.UIAxes_Viewer;
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

            id_z = this.nuid{zidx};
            if ~isempty(id_z)
                h_obj = cell(numel(id_z),1);
                % draw circle at the plane
                for cir_index = 1:numel(id_z)
                    cir_id = id_z(cir_index);

                    % skip the removed object
                    if cir_id == this.REMOVED_OBJ, continue; end

                    if ~isnan(cir_id)
                        self = images.roi.Circle(this.ax_caller,...
                            "Center",this.center{zidx}(cir_index,:),...
                            "Radius",this.radius{zidx}(cir_index), ...
                            "Color",this.color{zidx}(cir_index,:), ...
                            "LabelTextColor", this.VALID_FONT_COLOR, ...
                            "LineWidth",2,"MarkerSize",0.1,"Label",this.nlbl{zidx}(cir_index),...
                            "LabelVisible","on","FaceAlpha",0,"LabelAlpha",0,...
                            "Deletable",true,"UserData",[cir_index, zidx]);

                        % binding listener(will be auto removed if the roi was deleted)
                        addlistener(self, ...
                            'ROIMoved', @(src,~)this.roi_moved(src,[],self));

                        % binding menu: 修改ID
                        uimenu(self.ContextMenu, ...
                            "Text","修改ID",...
                            "MenuSelectedFcn",@(~,~)this.roi_modify_id([],[],self), ...
                            "Tag","IDROIContextMenuModify");

                        % binding menu: 修改标签
                        uimenu(self.ContextMenu, ...
                            "Text","修改标签",...
                            "MenuSelectedFcn",@(~,~)this.roi_modify_label([],[],self), ...
                            "Tag","LBLROIContextMenuModify");

                        % changed the default context menu: delete
                        self.ContextMenu.Children(end).MenuSelectedFcn ...
                            = @(~,~)this.roi_delete([],[],self);

                        h_obj{cir_index} = self;
                    end
                end
            else
                h_obj = {};
            end

            this.hobj = h_obj;
            status = this.STATUS_SUCCESS;
        end

        function status = Clear(this)
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
            this.nlbl{zidx}(cir_index+1,1) = string(cir_id);
            this.center{zidx}(cir_index+1, :) = cir_info(1:2);
            this.radius{zidx}(cir_index+1,1) = cir_info(3);

            this.rois = GenROI(this.nuid, this.center, this.radius);

            this.gen_colormap();

            status = this.STATUS_SUCCESS;
        end

        function lbl = getLabel(this, id)
            loc = cellfun(@(x)any(x==id), this.nuid, "UniformOutput",true);
            subid = this.nuid(loc);
            zidx = find(loc);
            cidx = cellfun(@(x)find(x==id), subid, "UniformOutput",true);

            lbl = this.nlbl{zidx(1)}(cidx(1));
        end

        function status = Refresh(this)
            % this function reorder objects identities
            % access each object and modify the identities
            
            % apply the refreshed identity
            [this.nuid, this.center, this.radius] ...
                = this.nuclears.ApplyRefreshTo(this.nuid, ...
                this.center, ...
                this.radius, ...
                true);

            % cover old names
            this.nlbl = cellfun(@(x)string(x),this.nuid,"UniformOutput",false);

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

        function status = recall_object(this, zidx, cir_index, cir_id_old, cir_id_new, cir_lbl_new)
            % this function rename an object in data already
            arguments
                this
                zidx        (1,1) double {mustBePositive, mustBeInteger}
                cir_index   (1,1) double  {mustBeNonNan} % empty for not needed: cir_id_old is not nan
                cir_id_old  (1,1) double
                cir_id_new  (1,1) double
                cir_lbl_new (1,1) string = ""
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
                        this.modify_label(cir_id_new, cir_lbl_new);
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
                this.modify_label(cir_id_new, cir_lbl_new);
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
                self.Label =  this.nlbl{zidx}(cidx);
                self.Color = this.color{zidx}(cidx, :);
            else
                warning("NuclearCtrlBot:invalidIdentity", ...
                    "Identity must be a positive integer.");
            end
        end

        function roi_modify_label(this, ~, ~, self)
            % generate a message box for new label acquirement
            cidx = self.UserData(1);
            zidx = self.UserData(2);    % [circle_index, z_index]
            id_cur = this.nuid{zidx}(cidx);
            lbl_cur = this.nlbl{zidx}(cidx);

            prompt = {'Enter new label:'};
            dlgtitle = 'remark';
            fieldsize = [1 45];
            definput = {lbl_cur};

            answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
            if isempty(answer), return; end

            lbl_new = string(answer{1});

            this.recall_object(zidx, cidx, id_cur, id_cur, lbl_new);

            self.Label = lbl_new;
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

        function modify_label(this, id, lbl)
            % find all identity equals to unid{zidx}(cidx), change the ulbl
            % as lbl, otherwise lbl is ""
            if lbl == ""
                return;
            else
                % modify all labels
                loc = cellfun(@(x)(id==x), this.nuid, "UniformOutput",false);
                for n = 1:numel(loc)
                    this.nlbl{n}(loc{n}) = lbl;
                end
            end
        end
    end
end