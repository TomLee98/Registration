classdef NuclearSplitter < handle
    %NUCLEARSPLITTER This class is nuclear splitter definition
    
    properties(Constant, Hidden)
        P_VALUE = 0.05;
        SHRINK_COEFFICIENT = 0.85;
        MIN_CIRCLE_RADIUS_PLANE = 4;
        MAX_SIZE_ALLDIM = 128;
    end

    properties(Access=private)
        %% input data
        segopts         % segment options
        caller          % the caller, must with SetProgressBar method

        %% output data
        center
        radius
        nid
        nuclears        % NuclearGroup object
        volume
    end

    properties(GetAccess=public, Dependent)
        Parent
        Center
        Radius
        Nid
    end

    methods
        function r = get.Parent(this)
            r = this.caller;
        end

        function r = get.Center(this)
            r = this.center;
        end

        function r = get.Radius(this)
            r = this.radius;
        end

        function r = get.Nid(this)
            r = this.nid;
        end
    end

    methods(Access=public)
        function this = NuclearSplitter(seg_opts, caller)
            %COUNTER A constructor
            arguments
                seg_opts  (1,1)     segopt
                caller    (1,1)     Register
            end

            this.segopts = seg_opts;
            this.caller = caller;
        end

        function Split(this, vol, volopts)
            % call this function for nuclear counting
            % do calculation on vol-volopts dataset
            switch this.segopts.method
                case "morphology"
                    [center_,radius_,id_] = split_morphology(this, vol, volopts);
                case "cellpose"
                    [center_,radius_,id_] = split_cellpose(this, vol, volopts);
                otherwise
            end
            
            this.center = center_;
            this.radius = radius_;
            this.nid = id_;
            this.volume = vol;
            this.nuclears = NuclearGroup(id_);
        end
    end

    methods(Access=private)
        function [center,radius,id] = split_morphology(this, vol, volopts)
            % This function uses morphology method to extract nuclears
            % preprocess the volume
            this.caller.SetProgressBar(0);

            bdbox = this.segopts.bdbox;
            % use watershed algorithm for objects segmentation
            bdbvol = cell(1, size(bdbox,1));
            for k = 1:numel(bdbvol)
                vol_k = imcrop3(vol, bdbox(k,:));     % crop volume

                % binarize the volume and resize to (128,?,?), where 128 is
                % maximum size of all dimension
                [bdbvol{k}, vopts] = NuclearSplitter.binarize(vol_k, volopts, this.segopts);

                this.caller.SetProgressBar(0.1*k/numel(bdbvol));
            end

            % utility options
            uopts = struct("pval",       this.P_VALUE, ...
                           "shrinkcoeff",this.SHRINK_COEFFICIENT);

            p = uopts.pval;

            r_range = [round((this.segopts.r_u+norminv(p/2)*this.segopts.r_s)/vopts.vr),...
                round((this.segopts.r_u-norminv(p/2)*this.segopts.r_s)/vopts.vr)];

            % parallel watershed on splitted volumes
            Ls = cell(1, numel(bdbvol));
            parfor k = 1:numel(bdbvol)
                % ============ Watershed for Initial Segmentation =============
                D = bwdist(~bdbvol{k}, "quasi-euclidean");
                D = -D;
                D(~bdbvol{k}) = Inf;
                Ls{k} = watershed(D, 26);   % rather than watershed_old
                Ls{k}(~bdbvol{k}) = 0;
                Ls{k} = imresize3(Ls{k}, bdbox(k,[5,4,6])+1, "Method","nearest"); %#ok<PFBNS>
            end

            % combine labels in L
            % assume that bdbox count is less than 256
            L = zeros(size(vol), "uint16");
            n_bdobj = 0;
            for k = 1:numel(bdbvol)
                % generate full L
                vbd = [bdbox(k,2), bdbox(k,2)+bdbox(k,5); ...  % rows
                       bdbox(k,1), bdbox(k,1)+bdbox(k,4); ...  % cols
                       bdbox(k,3), bdbox(k,3)+bdbox(k,6)];     % depth
                dnobj = max(Ls{k},[],"all");
                Ls{k}(Ls{k}~=0) = Ls{k}(Ls{k}~=0) + n_bdobj;
                L(vbd(1,1):vbd(1,2), vbd(2,1):vbd(2,2), vbd(3,1):vbd(3,2)) = Ls{k};
                n_bdobj = n_bdobj + dnobj;
            end
            L = imresize3(L, "Scale", [volopts.xRes, volopts.yRes, volopts.zRes]/vopts.vr, ...
                "Method", "nearest");

            % =============== Refine Objects by Criterion ================
            DIAMETER_MIN = 2*r_range(1);   % pixels
            DIAMETER_MAX = 2*r_range(2);

            stats = regionprops3(L, "Volume", "EquivDiameter", ...
                "Image", "SubarrayIdx","SurfaceArea");

            this.caller.SetProgressBar(0.5);

            for k = 1:size(stats, 1)
                sph = (4*pi)^(1/3)*(3*stats.Volume(k))^(2/3)/stats.SurfaceArea(k);

                % validate diameter and sphericity
                if stats.EquivDiameter(k) < DIAMETER_MIN ...
                        || stats.EquivDiameter(k) > DIAMETER_MAX ...
                        || sph < this.segopts.sph_th

                    % set as background
                    L(L==k) = 0;
                end

                this.caller.SetProgressBar(0.5+0.15*k/size(stats, 1));
            end

            % ============= Limit Objects Number by Maximum N ===========
            stats = regionprops3(L, "Volume");
            lbl_old = 1:numel(stats.Volume);
            lbl_old(stats.Volume==0) = [];
            stats(stats.Volume==0, :) = [];

            if numel(stats.Volume) > this.segopts.n_max
                dn = numel(stats.Volume) - this.segopts.n_max;
                % remove outliers
                switch this.segopts.rm_outlier
                    case "both"
                        nv = prctile(stats.Volume, [50*dn/numel(stats.Volume), ...
                            100-50*dn/numel(stats.Volume)]);
                        lbl_rm = lbl_old(stats.Volume<=nv(1)|stats.Volume>=nv(2));
                    case "left"
                        nv = prctile(stats.Volume, 100*dn/numel(stats.Volume));
                        lbl_rm = lbl_old(stats.Volume<=nv);
                    case "right"
                        nv = prctile(stats.Volume, 100-100*dn/numel(stats.Volume));
                        lbl_rm = lbl_old(stats.Volume>=nv);
                    otherwise
                end

                % remove labels
                for lbl_k = lbl_rm; L(L==lbl_k) = 0; end
            end

            this.caller.SetProgressBar(0.75);

            % ============ Reshape Objects for extraction ==============
            % re-sample (most down sampling as raw)
            L = imresize3(L, size(vol), "Method","nearest");

            % remap labels
            oldlbl = unique(L);
            newlbl = 0:numel(oldlbl)-1;
            for k = 1:numel(newlbl)
                L(L==oldlbl(k)) = newlbl(k);
            end

            % ========= Transform Label Matrix as Inner Storage ========
            % each c element: [x, y, warn_flag]
            center = cell(size(L, 3), 1);     % each element indicate circles in plane: center
            radius = cell(size(L, 3), 1);     % each element indicate circles in plane: radius
            id = cell(size(L, 3), 1);    % each element indicate circles in plane: identity

            se = strel("disk", ...
                max(1, round(max(vopts.xr_old, vopts.yr_old)/vopts.vr)));
            for n = 1:numel(id)
                % get regionprop at each plane and allocate identity
                Ln = L(:,:,n);
                Ln = imerode(Ln, se);   % imerode to split connected region because down sampling
                Vn = vol(:,:,n);

                % relabel stay order for regionprops linearity performance
                rLn = Ln;
                oldlbl = unique(Ln);
                newlbl = min(oldlbl):min(oldlbl)+numel(unique(Ln))-1;
                for k = 1:numel(newlbl)
                    rLn(Ln==oldlbl(k)) = newlbl(k);
                end

                stats = regionprops("table", rLn, Vn, "Centroid", ...
                    "EquivDiameter","Circularity","PixelIdxList");
                if ~isempty(stats)
                    % remove the lost label area or too small area
                    stats(stats.EquivDiameter<=2*this.MIN_CIRCLE_RADIUS_PLANE, :) = [];

                    center{n} = stats.Centroid;
                    radius{n} = stats.EquivDiameter/2*uopts.shrinkcoeff;
                    uncertain_id = (stats.Circularity < this.segopts.sph_th);

                    % get the center index of a region
                    idx = cellfun(@(x)x(round(numel(x)/2)), stats.PixelIdxList, ...
                        "UniformOutput",true);
                    id{n} = double(reshape(Ln(idx), [], 1));
                    id{n}(uncertain_id) = nan;

                    % remove conflict marker
                    id{n}(id{n}==0) = [];
                else
                    center{n} = double.empty(0,2);
                    radius{n} = double.empty(0,1);
                    id{n} = double.empty(0,1);
                end

                this.caller.SetProgressBar(0.75+0.15*n/numel(id));
            end

            % ================ Rearrange the Identities =================
            oldlbl = unique(cell2mat(id));
            oldlbl(isnan(oldlbl)) = [];
            newlbl = 1:numel(oldlbl);
            for k = 1:numel(id)
                if ~isempty(id{k})
                    [~, loc] = ismember(id{k}, oldlbl);
                    if any(~isnan(id{k}))
                        loc(loc==0) = [];
                        id{k}(~isnan(id{k})) = newlbl(loc);
                    end
                end

                this.caller.SetProgressBar(0.9+0.1*k/numel(id));
            end
        end

        function [center, radius, id] = split_cellpose(this, vol, volopts)
            % This function uses cellpose (machine learning) method to extract nuclears
            %% pre-process
            this.caller.SetProgressBar(0);
            bdvol = cell(1, size(this.segopts.bdbox,1));
            for k = 1:size(this.segopts.bdbox, 1)
                bdvol{k} = imcrop3(vol, this.segopts.bdbox(k,:));     % crop volume
            end

            slices = size(bdvol{1}, 3);

            for k = 1:numel(bdvol)
                % padding to square
                [height, width] = size(bdvol{k}, [1,2]);
                wp = max(height, width);
                if height == width
                    % already square image
                    dw(k,:) = [0, 0]; %#ok<AGROW>
                    dh(k,:) = [0, 0]; %#ok<AGROW>
                else
                    dw(k,:) = [floor((wp - width)/2), ceil((wp - width)/2)]; %#ok<AGROW>
                    dh(k,:) = [floor((wp - height)/2), ceil((wp - height)/2)]; %#ok<AGROW>
                end
                bdvol{k} = padarray(bdvol{k}, [dh(k,1), dw(k,1)], "replicate","pre");
                bdvol{k} = padarray(bdvol{k}, [dh(k,2), dw(k,2)], "replicate","post");

                % padding some slices for cubic
                rvol = uint8(rescale(bdvol{k},0,255));
                T = graythresh(rvol);
                mc_zv = reshape(sum(rvol-T*255,[1,2]),[],1);
                mc_z(k) = round(sum(mc_zv.*(0.5:1:slices-0.5)')./sum(mc_zv)); %#ok<AGROW>
                z_bdsiz = round(wp/volopts.zRes*volopts.xRes);
                dss = z_bdsiz - slices;
                ds(k,:) = [floor(dss/2), ceil(dss/2)]; %#ok<AGROW>
                if any(ds(k,:) < 0)
                    % remove some slices
                    rds(k,:) = [floor((z_bdsiz-1)/2), ceil((z_bdsiz-1)/2)]; %#ok<AGROW>
                    bdvol{k} = bdvol{k}(:,:,max(mc_z(k)-rds(k,1),1):min(mc_z(k)+rds(k,2),slices));
                else
                    bdvol{k} = padarray(bdvol{k}, [0, 0, ds(k,1)], "replicate","pre");
                    bdvol{k} = padarray(bdvol{k}, [0, 0, ds(k,2)], "replicate","post");
                end
                ps(k,:) = size(bdvol{k}); %#ok<AGROW>

                % resize to cuboid
                bdvol{k} = imresize3(bdvol{k}, wp*ones(1,3), "Method","linear");
                this.caller.SetProgressBar(0.1/numel(bdvol));
            end

            %% segmentation: call cellpose toolbox
            model = this.segopts.model;
            accalg = this.segopts.acceleration;
            cell_rk_th = this.segopts.cell_th;
            batch = this.segopts.batch_size;
            r_mean = this.segopts.r_mean/volopts.xRes;
            r_min = this.segopts.r_min/volopts.xRes;
            enh = this.segopts.enhance_ml;
            stage = this.segopts.stage;

            % run cellpose on CPU or GPU
            if gpuDeviceCount > 0
                cellpose_on_gpu();
            else
                cellpose_on_cpu();
            end

            for k = 1:numel(bdvol)
                % mass center reconstruction for better nuclear segmentation
                % foreground, and decrease the effect of residual motion
                rp = regionprops3(labels{k}, "Centroid");
                if ~isempty(rp)
                    mask = false(size(bdvol{k}));
                    cc = round(rp.Centroid(:,[2,1,3]));
                    mask(sub2ind(size(bdvol{k}), cc(:,1), cc(:,2), cc(:,3))) = true;

                    % dilate for nuclear reconstruction
                    mask = imdilate(mask, strel("sphere", ...
                        max(1, round(this.segopts.r_mean/volopts.xRes))));

                    % refine labels
                    labels{k} = labels{k}.*mask; %#ok<AGROW>
                end
            end

            this.caller.SetProgressBar(0.85);

            %% post-process
            n_bdobj = 0;
            % extract objects center, radius and identity
            center = cell(slices, 1);
            radius = cell(slices, 1);
            id = cell(slices, 1);

            for k = 1:numel(bdvol)
                % recovery volume as raw
                labels{k} = imresize3(labels{k}, ps(k,:), "Method","nearest");

                if any(ds(k,:) < 0)
                    % padding slices
                    zb = max(mc_z(k)-rds(k,1),1); ze = min(mc_z(k)+rds(k,2),slices);
                    labels{k} = padarray(labels{k}, [0,0,zb-1], "replicate", "pre");
                    labels{k} = padarray(labels{k}, [0,0,slices-ze], "replicate", "post");
                    labels{k} = labels{k}(dh(k,1)+1:end-dh(k,2), dw(k,1)+1:end-dw(k,2), :);
                else
                    % remove padding slices
                    labels{k} = labels{k}(dh(k,1)+1:end-dh(k,2), dw(k,1)+1:end-dw(k,2), ...
                        ds(k,1)+1:end-ds(k,2));
                end

                % remap labels after volume resizing
                old_label = unique(labels{k}); old_label(old_label==0) = [];
                for m = 1:numel(old_label)
                    labels{k}(labels{k}==old_label(m)) = m;
                end

                % extract objects on each plane
                for s = 1:slices
                    rp = regionprops("table", labels{k}(:,:,s), "Centroid", "EquivDiameter");
                    if ~isempty(rp)
                        id{s} = [id{s}; find(rp.EquivDiameter>eps) + n_bdobj];
                        center{s} = [center{s}; rp.Centroid(rp.EquivDiameter>eps,:) ...
                            + this.segopts.bdbox(k,1:2)];
                        radius{s} = [radius{s}; rp.EquivDiameter(rp.EquivDiameter>eps)/2];
                    else
                        id{s} = [id{s}; double.empty(0,1)];
                        center{s} = [center{s}; double.empty(0,2)];
                        radius{s} = [radius{s}; double.empty(0,1)];
                    end
                    this.caller.SetProgressBar(0.85+0.15*s*k/slices/numel(bdvol));
                end

                n_bdobj = n_bdobj + max(labels{k},[],"all");   % update objects number
            end

            % parallel on CPU
            function cellpose_on_cpu()
                % for loop is faster because pyenv takes over CPU parallel
                for n = 1:numel(bdvol)
                    switch model
                        case {'cyto', 'CP', 'nuclei'}
                            cp = cellpose("Model", model, ...
                                "Acceleration", accalg);
                            labels{n} = segmentCells3D(cp, bdvol{n}, ...
                                "CellThreshold", cell_rk_th, ...
                                "GPUBatchSize", batch, ...
                                "ImageCellDiameter", round(2*r_mean), ...
                                "MinVolume",round(4/3*pi*r_min.^3), ...
                                "TileAndAugment", enh);
                        case "KC"
                            [mpath, ~, ~] = fileparts(mfilename("fullpath"));
                            mpath = [mpath, filesep, 'custom_models']; %#ok<AGROW>
                            cp = cellpose("Model", mpath + filesep + model + "_" + stage, ...
                                "Acceleration", accalg);
                            labels{n} = segmentCells3D(cp, bdvol{n}, ...
                                "CellThreshold", cell_rk_th, ...
                                "GPUBatchSize", batch, ...
                                "ImageCellDiameter", round(2*r_mean), ...
                                "MinVolume",round(4/3*pi*r_min.^3), ...
                                "TileAndAugment", enh);
                        otherwise
                    end

                    terminate(pyenv);

                    this.caller.SetProgressBar(0.1+n/numel(bdvol)*0.7);
                end
            end

            % parallel on GPU
            function cellpose_on_gpu()
                % foreground client generate job
                parfor n = 1:numel(bdvol)
                    job_n = parfeval(@NuclearSplitter.cellpose_on_each_gpu, 1, ...
                        model, bdvol{n}, cell_rk_th, stage, batch, r_mean, r_min, enh);

                    % get background clients results
                    labels{n} = fetchOutputs(job_n);
                end

                this.caller.SetProgressBar(0.8);
            end
        end
    end

    methods (Static)
        % rehist: re histogram the bright to enhance contrast
        function [bvol, opts] = binarize(vol, volopts, segopts)
            opts.xr_old = volopts.xRes;
            opts.yr_old = volopts.yRes;
            opts.zr_old = volopts.zRes;
            opts.w = volopts.width;
            opts.h = volopts.height;
            opts.s = volopts.slices;
            rs = [volopts.xRes, volopts.yRes, volopts.zRes];
            opts.vr = min(rs);
            opts.vr = opts.vr*max(rs/opts.vr.*size(vol))/NuclearSplitter.MAX_SIZE_ALLDIM;

            % remap to 0-255 as uint8
            H_MIN = min(vol, [], "all");
            vol = vol - H_MIN;
            pvol = uint8(rescale(single(vol).^segopts.gamma, 0, 255));

            % use filter for smooth volume
            pvol = medfilt3(pvol, [3,3,1], "replicate");
            pvol = imgaussfilt3(pvol, [1,1,1]);

            % up-sampling for better watershed processing
            pvol = imresize3(pvol, "Scale", ...
                [volopts.xRes, volopts.yRes, volopts.zRes]/opts.vr, ...
                "Method", "linear");

            % binarize foreground
            if segopts.auto_fg
                bvol = imbinarize(pvol, "global");
            else
                bvol = (pvol > segopts.fg_th);
            end
        end

        function label = cellpose_on_each_gpu(model, bdvol, cth, stage, batch, r_mean, r_min, enh)
            switch model
                case {'cyto', 'CP', 'nuclei'}
                    cp = cellpose("Model", model);
                    label = segmentCells3D(cp, bdvol, ...
                        "CellThreshold", cth, ...
                        "GPUBatchSize", batch, ...
                        "ImageCellDiameter", round(2*r_mean), ...
                        "MinVolume",round(4/3*pi*r_min.^3), ...
                        "TileAndAugment", enh);
                case "KC"
                    [mpath, ~, ~] = fileparts(mfilename("fullpath"));
                    mpath = [mpath, filesep, 'custom_models'];
                    cp = cellpose("Model", mpath + filesep + model + "_" + stage);
                    label = segmentCells3D(cp, bdvol, ...
                        "CellThreshold", cth, ...
                        "GPUBatchSize", batch, ...
                        "ImageCellDiameter", round(2*r_mean), ...
                        "MinVolume",round(4/3*pi*r_min.^3), ...
                        "TileAndAugment", enh);
                otherwise
            end

            terminate(pyenv);
        end
    end
end