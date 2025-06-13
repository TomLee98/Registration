classdef regrspt < handle
    %REGRSPTS This class defines Restore-point object, which restored the
    % reg3D processed data

    properties (Constant, Hidden)
        RESTORE_PARALLEL_POOL_SIZE = 24
        VALID_FIELD_NAME = ["text_time", "text_meta", "fixdef", "mask_reg"]

    end
    
    properties (GetAccess = private, SetAccess = immutable)
        args    (1,1)   struct                                                              % argument struct, with field {"REGISTER", "SEGMENT", "CROP", "LOAD"}  
                                                                                            % value could be regopt object, segopt object, struct(for crop) or []
        cdim    (1,1)   string  {mustBeMember(cdim, ["XY","Z","T",""])}                     % crop dimension 
        dptr    (1,1)   regmov                                           = regmov.empty()   % regmov object as stored data
        optr    (1,1)   string                                                              % operator, could be "CROP"/"REGISTER"/"SEGMENT"
        others  (1,1)   struct                                                              % other components binding on reg3d
        tfs     (:,1)   cell                                                                % additional parameter as transformation objects, if optr as "REGISTER"
        ncbptr  (1,:)                                                                       % NuclearCtrlBot object
    end

    properties (Access = private)

    end

    properties (Access = public, Dependent)
        Arguments       % ___/get, 1-by-1 struct with regopt, segopt, struct and []
        CropDim         % ___/get, 1-by-1 string indicate crop dimension
        CellsCount      % ___/get, 1-by-1 nonnegative integer indicates cells count
        File            % ___/get, 1-by-1 string indicate inner data file
        ImageDim        % ___/get, 1-by-5 positive integer as [X,Y,C,Z,T] dimensions
        IsDistributed   % set/get, 1-by-1 logical, indicate the storage data location
        Operator        % ___/get, 1-by-1 string, operator indicator
        Others          % ___/get, 1-by-1 struct with binding data
        Segmentor       % ___/get, 1-by-1 NuclearCtrlBot or []
    end

    methods

        function r = get.Arguments(this)
            r = this.args;
        end

        function r = get.CropDim(this)
            if this.cdim ~= ""
                r = this.cdim;
            else
                r = string([]); % not applicable
            end
        end

        function r = get.File(this)
            if isempty(this.dptr)
                r = struct("File","", "Size",nan);  % not applicable
            else
                r = struct("File", this.dptr.Location, ...
                           "Size", this.dptr.Bytes);
            end
        end

        function r = get.CellsCount(this)
            if isempty(this.ncbptr)
                r = nan;        % not applicable
            else
                r = this.ncbptr.MaxN;
            end
        end

        function r = get.ImageDim(this)
            md = this.dptr.MetaData;
            r = [md.width, md.height, md.channels, md.slices, md.frames];
        end

        function r = get.IsDistributed(this)
            r = isfile(this.dptr.Location);
        end

        function set.IsDistributed(this, r)
            arguments
                this
                r       (1,1)   logical
            end

            if r == true
                if ~isempty(this.dptr), this.dptr.spread(); end
            else
                if ~isempty(this.dptr), this.dptr.gather(); end
            end
        end

        function r = get.Operator(this)
            r = this.optr;
        end

        function r = get.Others(this)
            r = this.others;
        end

        function r = get.Segmentor(this)
            r = this.ncbptr;
        end
    end
    
    methods (Access = public)
        function this = regrspt(optr, args, varargin)
            %REGRSPTS A Constructor
            p = inputParser();
            p.StructExpand = false;
            addRequired(p, "optr", @(x)validateattributes(x, "string", "scalar"));
            addRequired(p, "args", @(x)validateattributes(x, "struct", "scalar"));
            addParameter(p, "CropDim", "");
            addParameter(p, "Data", regmov.empty());
            addParameter(p, "Others", struct("text_time",[], "text_meta",[], "fixdef",[], "mask_reg",[]));
            addParameter(p, "Segmentor", []);
            addParameter(p, "Transform", {});

            p.parse(optr, args, varargin{:});

            if ~ismember(p.Results.optr, [constdef.OP_LOAD, constdef.OP_CROP, ...
                    constdef.OP_REGISTER, constdef.OP_SEGMENT])
                throw(MException("regrspt:invalidOperator", ...
                    "RestorePoint object only support <LOAD>, <CROP>, <REGISTER>, <SEGMENT>"));
            else
                this.optr = p.Results.optr;
                this.args = p.Results.args;
            end

            this.cdim = p.Results.CropDim;
            this.dptr = p.Results.Data;
            
            switch this.optr
                case constdef.OP_REGISTER
                    switch this.args.(constdef.OP_REGISTER).Mode
                        case "Global"
                            this.tfs = p.Results.Transform(:,1);
                        case "Local"
                            this.tfs = p.Results.Transform(:,2);
                        case "Manual"
                            this.tfs = p.Results.Transform(:,3);
                        otherwise
                    end
                otherwise
            end

            if ~isempty(setdiff(string(fieldnames(p.Results.Others)), this.VALID_FIELD_NAME))
                throw(MException("regrspt:invalidOthers", ...
                    "Invalid field in struct 'Others'."));
            end

            this.others = p.Results.Others;

            this.ncbptr = p.Results.Segmentor;
        end
        
        function dpost = getData(this, dpre)
            % getData This function get the restored data from this
            arguments
                this 
                dpre (1,1)  regmov = regmov.empty()
            end

            switch this.optr
                case constdef.OP_LOAD
                    dpost = this.restore_load();
                case constdef.OP_CROP
                    dpost = this.restore_crop(dpre);
                case constdef.OP_REGISTER
                    dpost = this.restore_register(dpre);
                case constdef.OP_SEGMENT
                    dpost = this.restore_segment(dpre);
                otherwise
                    throw(MException("regrspt:invalidOperator", ...
                        "RestorePoint object only support <CROP>, <REGISTER>, <SEGMENT>"));
            end
        end

        % invoke this to free related data
        function delete(this)
            % ~
            % set the dptr "Retain-Cache" to false
            this.dptr.RetainCache = false;  % 
            delete(this.dptr);

            delete(this.ncbptr);
        end
    end

    methods (Access = private)
        function dpost = restore_load(this)
            if isempty(this.dptr)
                throw(MException("regrspt:invalidInvoke", ...
                    "No stored data can be restored."));
            else
                dpost = this.dptr;      % pass
            end
        end

        function dpost = restore_crop(this, dpre)
            if isempty(dpre)            % no given data
                dpost = get_storage_data(this);
            else
                % do crop on dpre
                switch this.cdim
                    case "XY"
                        dpost = dpre.vcrop("XY", this.args.(constdef.OP_CROP).xy');
                    case "Z"
                        dpost = dpre.vcrop("Z", this.args.(constdef.OP_CROP).z);
                    case "T"
                        dpost = dpre.tcrop(this.args.(constdef.OP_CROP).f);
                    otherwise
                        throw(MException("regrspt:invalidCrop", ...
                            "No valid crop dimension."));
                end
            end
        end

        function dpost = restore_register(this, dpre)
            if isempty(dpre)
                dpost = get_storage_data(this);
            else
                % do imwarp on dpre
                dpost = dpre.copy();    % deep copy
                rref = imref3d([dpre.MetaData.Height, dpre.MetaData.Width, dpre.MetaData.Slices], ...
                    movsrc.MetaData.xRes, movsrc.MetaData.yRes, movsrc.MetaData.zRes);
                itpalg = this.args.(constdef.OP_REGISTER).Interp;
                frs = ceil(1:this.RESTORE_PARALLEL_POOL_SIZE:dpre.MetaData.Frames+1);
                try
                    delete(gcp("nocreate"));
                    pp = parpool("Threads", this.RESTORE_PARALLEL_POOL_SIZE);
                catch
                    throw(MException("regrspt:invalidParallelToolboxState", ...
                        "Can't reset parpool or change the workers number."));
                end
                
                % blocked data loading 
                for r = 1:numel(frs)-1
                    % load data
                    regfrs = frs(r):frs(r+1)-1;
                    fmode = ["none", string(regfrs).join(",")];
                    avol_sc = grv(dpre, fmode, this.args.(constdef.OP_REGISTER).SC);
                    avol_fc = grv(dpre, fmode, this.args.(constdef.OP_REGISTER).FC);
                    tfss = this.tfs(regfrs);

                    switch this.args.(constdef.OP_REGISTER).Algorithm
                        case "TCREG"
                            switch this.args.(constdef.OP_REGISTER).Mode
                                case "Global"
                                    parfor m = 1:this.RESTORE_PARALLEL_POOL_SIZE
                                        avol_sc_m = avol_sc(:,:,:,m);
                                        avol_fc_m = avol_fc(:,:,:,m);
                                        fival_sc = mean(avol_sc_m(:,[1,end],:),"all");
                                        fival_fc =  mean(avol_fc_m(:,[1,end],:),"all");
                                        avol_sc(:,:,:,m) = imwarp(avol_sc_m, ...
                                            tfss{m}, itpalg,  "OutputView", rref, ...
                                            "FillValues",fival_sc);
                                        avol_fc(:,:,:,m) = imwarp(avol_fc_m, ...
                                            tfss{m}, itpalg, "OutputView", rref, ...
                                            "FillValues",fival_fc);
                                    end
                                case "Local"
                                    parfor m = 1:this.RESTORE_PARALLEL_POOL_SIZE
                                        avol_sc_m = avol_sc(:,:,:,m);
                                        avol_fc_m = avol_fc(:,:,:,m);
                                        fival_sc = mean(avol_sc_m(:,[1,end],:),"all");
                                        fival_fc =  mean(avol_fc_m(:,[1,end],:),"all");
                                        avol_sc(:,:,:,m) = imwarp(avol_sc_m, ...
                                            tfss{m}, itpalg, "FillValues",fival_sc);
                                        avol_fc(:,:,:,m) = imwarp(avol_fc_m, ...
                                            tfss{m}, itpalg, "FillValues",fival_fc);
                                    end
                                otherwise
                            end

                            % save warped data
                            srv(dpost, avol_sc, fmode, this.args.(constdef.OP_REGISTER).SC);
                            srv(dpost, avol_fc, fmode, this.args.(constdef.OP_REGISTER).FC);
                        case "OCREG"
                            switch this.args.(constdef.OP_REGISTER).Mode
                                case "Global"
                                    parfor m = 1:this.RESTORE_PARALLEL_POOL_SIZE
                                        avol_fc_m = avol_fc(:,:,:,m);
                                        fival_fc =  mean(avol_fc_m(:,[1,end],:),"all");
                                        avol_fc(:,:,:,m) = imwarp(avol_fc_m, ...
                                            rref, tfss{m}, itpalg, "OutputView", rref, ...
                                            "FillValues",fival_fc);
                                    end
                                otherwise
                            end

                            % save warped data
                            srv(dpost, avol_fc, fmode, this.args.(constdef.OP_REGISTER).FC);

                        case "MANREG"
                            % must be dual-color imaging
                            parfor m = 1:this.RESTORE_PARALLEL_POOL_SIZE
                                avol_sc_m = avol_sc(:,:,:,m);
                                avol_fc_m = avol_fc(:,:,:,m);
                                fival_sc = mean(avol_sc_m(:,[1,end],:),"all");
                                fival_fc =  mean(avol_fc_m(:,[1,end],:),"all");
                                avol_sc(:,:,:,m) = imwarp(avol_sc_m, ...
                                    tfss{m}, itpalg,  "OutputView", rref, ...
                                    "FillValues",fival_sc);
                                avol_fc(:,:,:,m) = imwarp(avol_fc_m, ...
                                    tfss{m}, itpalg, "OutputView", rref, ...
                                    "FillValues",fival_fc);
                            end

                            % save warped data
                            srv(dpost, avol_sc, fmode, this.args.(constdef.OP_REGISTER).SC);
                            srv(dpost, avol_fc, fmode, this.args.(constdef.OP_REGISTER).FC);
                        otherwise
                    end
                end

                delete(pp);

                switch this.args.(constdef.OP_REGISTER).Algorithm
                    case "TCREG"
                        switch this.args.(constdef.OP_REGISTER).Mode
                            case "Global"
                                dpost.Transformation(:,1) = this.tfs;
                            case "Local"
                                dpost.Transformation(:,2) = this.tfs;
                        end
                    case "OCREG"
                        dpost.Transformation(:,1) = this.tfs;
                    case "MANREG"
                        dpost.Transformation(:,3) = this.tfs;
                    otherwise
                end
            end
        end

        function dpost = restore_segment(this, dpre) %#ok<INUSD>
            % pass
            dpost = dpre;
        end

        function r = get_storage_data(this)
            if isempty(this.dptr)
                throw(MException("regrspt:invalidInvoke", ...
                    "No stored data can be restored."));
            else
                r = this.dptr;      % just copy
            end
        end
    end

end
