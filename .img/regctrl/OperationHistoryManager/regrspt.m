classdef regrspt < handle
    %REGRSPTS This class defines Restore-point object, which restored the
    % reg3D processed data

    properties (Constant, Hidden)
        RESTORE_PARPOOL_SIZE = ceil(GetCPUWorkersMaxN([],[])/2);
        VALID_FIELD_NAME = ["text_time", "text_meta", "fixdef", "frames_reg", ...
            "mask_reg", "mode_zproj"]
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
            addParameter(p, "Others", struct("text_time",[], "text_meta",[], ...
                                             "fixdef",[], "frames_reg",[], ...
                                             "mask_reg",[], "mode_zproj",[]));
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
                        case "global"
                            this.tfs = p.Results.Transform(:,1);
                        case "local"
                            this.tfs = p.Results.Transform(:,2);
                        case "manual"
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
                        "RestorePoint object only support <LOAD>, <CROP>, <REGISTER>, <SEGMENT>."));
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

                pp = parpool("Threads", this.RESTORE_PARPOOL_SIZE);

                switch this.args.(constdef.OP_REGISTER).Algorithm
                    case "TCREG"
                        switch this.args.(constdef.OP_REGISTER).Mode
                            case "global"
                                recover_reg_tc_global(this);
                                dpost.Transformation(:,1) = this.tfs;
                            case "local"
                                recover_reg_tc_local(this);
                                dpost.Transformation(:,2) = this.tfs;
                            otherwise
                        end
                    case "OCREG"
                        switch this.args.(constdef.OP_REGISTER).Mode
                            case "global"
                                recover_reg_oc_global(this);
                                dpost.Transformation(:,1) = this.tfs;
                            otherwise
                        end
                    case "MANREG"
                        recover_reg_tc_manual(this);
                        dpost.Transformation(:,3) = this.tfs;
                    otherwise
                end

                delete(pp);

            end

            %% NESTED HELPER FUNCTIONS
            function recover_reg_tc_global(this)
                % This function recover tcreg global registration results
                rref = imref3d([dpre.MetaData.height, dpre.MetaData.width, dpre.MetaData.slices], ...
                    dpre.MetaData.xRes, dpre.MetaData.yRes, dpre.MetaData.zRes);
                itpalg = this.args.(constdef.OP_REGISTER).Options.Interp;

                fridx = str2num(this.others.frames_reg(1).replace("end", ...
                    string(dpre.MetaData.frames))); %#ok<ST2NM>
                fridxx = find(cellfun(@(x)~isempty(x), this.tfs, "UniformOutput",true));
                fridx = intersect(fridx, fridxx);

                RPS = this.RESTORE_PARPOOL_SIZE;
                if mod(numel(fridx), RPS) ~= 0
                    frs = [1:RPS:RPS*floor(numel(fridx)/RPS)+1, numel(fridx)+1];
                else
                    frs = 1:RPS:RPS*floor(numel(fridx)/RPS)+1;
                end

                for r = 1:numel(frs)-1
                    % load data
                    regfrs = fridx(frs(r):frs(r+1)-1);
                    fmode = ["none", string(regfrs).join(",")];
                    avol_sc = grv(dpre, fmode, this.args.(constdef.OP_REGISTER).Options.SC);
                    avol_fc = grv(dpre, fmode, this.args.(constdef.OP_REGISTER).Options.FC);
                    tfss = this.tfs(regfrs);

                    % imwarp
                    parfor m = 1:numel(regfrs)
                        avol_sc_m = avol_sc(:,:,:,m);
                        avol_fc_m = avol_fc(:,:,:,m);
                        fival_sc = mean(avol_sc_m(:,[1,end],:),"all");
                        fival_fc =  mean(avol_fc_m(:,[1,end],:),"all");
                        avol_sc(:,:,:,m) = imwarp(avol_sc_m, rref,...
                            tfss{m}, itpalg,  "OutputView", rref, ...
                            "FillValues",fival_sc);
                        avol_fc(:,:,:,m) = imwarp(avol_fc_m, rref,...
                            tfss{m}, itpalg, "OutputView", rref, ...
                            "FillValues",fival_fc);
                    end

                    % save warped data
                    srv(dpost, avol_sc, fmode, this.args.(constdef.OP_REGISTER).Options.SC);
                    srv(dpost, avol_fc, fmode, this.args.(constdef.OP_REGISTER).Options.FC);
                end
            end

            function recover_reg_tc_local(this)
                % This function recover tcreg local registration results
                itpalg = this.args.(constdef.OP_REGISTER).Options.Interp;

                fridx = str2num(this.others.frames_reg(1).replace("end", ...
                    string(dpre.MetaData.frames))); %#ok<ST2NM>
                fridxx = find(cellfun(@(x)~isempty(x), this.tfs, "UniformOutput",true));
                fridx = intersect(fridx, fridxx);
                RPS = this.RESTORE_PARPOOL_SIZE;
                if mod(numel(fridx), RPS) ~= 0
                    frs = [1:RPS:RPS*floor(numel(fridx)/RPS)+1, numel(fridx)];
                else
                    frs = 1:RPS:RPS*floor(numel(fridx)/RPS)+1;
                end

                for r = 1:numel(frs)-1
                    % load data
                    regfrs = fridx(frs(r):frs(r+1)-1);
                    fmode = ["none", string(regfrs).join(",")];
                    avol_sc = grv(dpre, fmode, this.args.(constdef.OP_REGISTER).Options.SC);
                    avol_fc = grv(dpre, fmode, this.args.(constdef.OP_REGISTER).Options.FC);
                    tfss = this.tfs(regfrs);

                    % imwarp
                    parfor m = 1:numel(regfrs)
                        if ~isempty(tfss{m})
                            avol_sc_m = avol_sc(:,:,:,m);
                            avol_fc_m = avol_fc(:,:,:,m);
                            fival_sc = mean(avol_sc_m(:,[1,end],:),"all");
                            fival_fc =  mean(avol_fc_m(:,[1,end],:),"all");
                            avol_sc(:,:,:,m) = imwarp(avol_sc_m, ...
                                tfss{m}, itpalg, "FillValues",fival_sc);
                            avol_fc(:,:,:,m) = imwarp(avol_fc_m, ...
                                tfss{m}, itpalg, "FillValues",fival_fc);
                        end
                    end

                    % save warped data
                    srv(dpost, avol_sc, fmode, this.args.(constdef.OP_REGISTER).Options.SC);
                    srv(dpost, avol_fc, fmode, this.args.(constdef.OP_REGISTER).Options.FC);
                end
            end

            function recover_reg_oc_global(this)
                % This function recover ocreg global registration results
                rref = imref3d([dpre.MetaData.height, dpre.MetaData.width, dpre.MetaData.slices], ...
                    dpre.MetaData.xRes, dpre.MetaData.yRes, dpre.MetaData.zRes);
                itpalg = this.args.(constdef.OP_REGISTER).Options.Interp;

                fridx = str2num(this.others.frames_reg(1).replace("end", ...
                    string(dpre.MetaData.frames))); %#ok<ST2NM>
                fridxx = find(cellfun(@(x)~isempty(x), this.tfs, "UniformOutput",true));
                fridx = intersect(fridx, fridxx);
                RPS = this.RESTORE_PARPOOL_SIZE;
                if mod(numel(fridx), RPS) ~= 0
                    frs = [1:RPS:RPS*floor(numel(fridx)/RPS)+1, numel(fridx)];
                else
                    frs = 1:RPS:RPS*floor(numel(fridx)/RPS)+1;
                end

                for r = 1:numel(frs)-1
                    % load data
                    regfrs = fridx(frs(r):frs(r+1)-1);
                    fmode = ["none", string(regfrs).join(",")];
                    avol_fc = grv(dpre, fmode, this.args.(constdef.OP_REGISTER).Options.FC);
                    tfss = this.tfs(regfrs);

                    % imwarp
                    parfor m = 1:numel(regfrs)
                        avol_fc_m = avol_fc(:,:,:,m);
                        fival_fc =  mean(avol_fc_m(:,[1,end],:),"all");
                        avol_fc(:,:,:,m) = imwarp(avol_fc_m, rref, ...
                            tfss{m}, itpalg, "OutputView", rref, ...
                            "FillValues",fival_fc);
                    end

                    % save warped data
                    srv(dpost, avol_fc, fmode, this.args.(constdef.OP_REGISTER).Options.FC);
                end
            end

            function recover_reg_tc_manual(this)
                % This function recover manreg global registration results
                rref = imref3d([dpre.MetaData.height, dpre.MetaData.width, dpre.MetaData.slices], ...
                    dpre.MetaData.xRes, dpre.MetaData.yRes, dpre.MetaData.zRes);
                itpalg = this.args.(constdef.OP_REGISTER).Options.Interp;

                % load data
                regfrs = str2double(this.others.frames_reg(2));
                fmode = ["none", string(regfrs).join(",")];
                avol_sc = grv(dpre, fmode, this.args.(constdef.OP_REGISTER).Options.SC);
                avol_fc = grv(dpre, fmode, this.args.(constdef.OP_REGISTER).Options.FC);
                tfss = this.tfs{regfrs};

                % imwarp
                fival_sc = mean(avol_sc(:,[1,end],:),"all");
                fival_fc =  mean(avol_fc(:,[1,end],:),"all");
                avol_sc = imwarp(avol_sc, rref, tfss, itpalg,  "OutputView", rref, ...
                    "FillValues",fival_sc);
                avol_fc = imwarp(avol_fc, rref, tfss, itpalg, "OutputView", rref, ...
                    "FillValues",fival_fc);

                % save warped data
                srv(dpost, avol_sc, fmode, this.args.(constdef.OP_REGISTER).Options.SC);
                srv(dpost, avol_fc, fmode, this.args.(constdef.OP_REGISTER).Options.FC);
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
