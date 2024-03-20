classdef regmask < handle
    %REGMASK This class is mask object defination, with basic operation on
    %mask

    properties(Constant, Hidden)
        REMOVE_THRESHOLD = 0.5
    end
    
    properties(Access=private, Hidden)
        bdbox           % s-by-6 matrix, each row is components bounding box, without background
        mask_vol        % m-by-n-by-p label volume
        mask_roi        % p+1-by-1 cell array, each item for one plane, first item with background
        mask_lbl        % 1-by-s string array, each item is a component label, without background
        mask_comps      % 1-by-s array, indicates the components label in mask_vol, without background
    end

    properties(GetAccess=public, Dependent)
        MaskVol         % get,      item is 3d labeled array
        MaskROI         % get,      item shape as: [cx, cy, r, id]
        MaskLabel       % get/set,  item is string
        MaskComps       % get,      the components identities array
        NumComps        % get,      number of components in mask
    end
    
    methods
        function this = regmask(mask_, label_)
            %REGMASK A constructor
            arguments
                mask_   {mustBeMask}
                label_  (1,:)   string  = ""
            end

            % transform mask and stored
            if iscell(mask_)
                this.mask_roi = mask_;
                if "" ~= label_
                    this.mask_vol = regmask.roi2vol(mask_);
                    this.mask_lbl = label_;
                else
                    [this.mask_vol, this.mask_lbl] = regmask.roi2vol(mask_);
                end
            else
                this.mask_vol = mask_;
                if "" ~= label_
                    this.mask_roi = regmask.vol2roi(mask_);
                    this.mask_lbl = label_;
                else
                    [this.mask_roi, this.mask_lbl] = regmask.vol2roi(mask_);
                end
            end

            %
            comps = unique(this.mask_vol, "sorted");
            comps(0 == comps) = [];
            this.mask_comps = reshape(comps, 1, []);

            %
            this.calc_bounding_box();
        end

        function r = get.MaskVol(this)
            r = this.mask_vol;
        end
        
        function r = get.MaskROI(this)
            r = this.mask_roi(2:end);   % skip the background
        end

        function r = get.MaskLabel(this)
            r = this.mask_lbl;
        end

        function set.MaskLabel(this, r_)
            arguments
                this
                r_      (1,:)   string
            end

            if numel(this.mask_lbl) ~= numel(r_)
                throw(MException("regmask:setMaskLabel", ...
                    "Invalid labels counts."));
            end

            this.mask_lbl = r_;
        end

        function r = get.MaskComps(this)
            r = this.mask_comps;
        end

        function r = get.NumComps(this)
            r = numel(this.mask_comps);
        end

        function r = size(this, dim)
            arguments
                this
                dim     (1,:)   double  {mustBeMember(dim, [1,2,3])} = [1,2,3]
            end
            r = size(this.mask_vol, dim);
        end

        function mask_obj = crop(this, type_, r_)
            % This function auto update components in mask when do crop on
            % the regmask object
            arguments
                this
                type_   (1,1)   string  {mustBeMember(type_, ["XY","Z"])}
                r_      (:,2)   {mustBePositive, mustBeInteger}
            end
            
            % if half size of an object is out of cropped boundary, then 
            % it will be remove from the new mask
            switch type_
                case "XY"
                    if size(r_, 1) ~= 2
                        throw(MException("regmask:crop:invalidRange", ...
                            "XY crop requires two dimension inputs."));
                    end

                    % crop the volume
                    this.mask_vol = this.mask_vol(r_(2,1):r_(2,2), r_(1,1):r_(1,2), :);

                    % loop at each plane, mark the cropped object
                    pcomp_rm = cell(numel(this.mask_roi)-1, 1);
                    for k = 1:numel(pcomp_rm)
                        if ~isempty(this.mask_roi{k+1})
                            bd = this.getBoundingBox(this.mask_roi{k+1}(:,4));
                            flag_kp = regmask.overlap_gt(r_, bd);
                            pcomp_rm{k} = this.mask_roi{k+1}(~flag_kp, 4);
                        end
                    end

                    % adjust comp_rm for objects to be removed
                    pcomp_rm = cell2mat(pcomp_rm);
                    comps_rm = unique(pcomp_rm,"sorted");
                    nz = histcounts(pcomp_rm, "BinMethod","integers");
                    nz(nz==0) = [];
                    nz = reshape(nz, [], 1);

                    bd = this.getBoundingBox(comps_rm);

                    % if number of lost slices less than threshold
                    % keep them stay, others to be removed
                    comps_rm(nz./bd(:,6) < this.REMOVE_THRESHOLD) = [];

                    comps_keep = setdiff(this.mask_comps, comps_rm);

                    mask_obj = this.extract(comps_keep);

                case "Z"
                    if size(r_, 1) ~= 1
                        throw(MException("regmask:crop:invalidRange", ...
                            "Z crop requires one dimension input."));
                    end

                    % crop the volume
                    this.mask_vol = this.mask_vol(:, :, r_(1):r_(2));

                    % skip first plane: indicate whole volume size
                    z_rm = setdiff(2:numel(this.mask_roi), r_(1)+1:r_(2)+1);
                    planes_rm =  this.mask_roi(z_rm);
                    pcomp_rm = cellfun(@(x)x(:, 4), planes_rm, "UniformOutput",false);
                    pcomp_rm = unique(cell2mat(pcomp_rm), "sorted");
                    
                    % adjust the objects to be removed
                    bdbox_rm = this.getBoundingBox(pcomp_rm);

                    comps_rm = [];
                    for k = 1:size(bdbox_rm, 1)
                        hb = round(bdbox_rm(k,3):bdbox_rm(k,3)+bdbox_rm(k,6)-1);
                        ratio_rm = ...
                            numel(intersect(hb, z_rm)) / bdbox_rm(k,6);

                        % exceed half size, remove it
                        if ratio_rm > this.REMOVE_THRESHOLD
                            % mark the removed object
                            comps_rm = [comps_rm, pcomp_rm(k)]; %#ok<AGROW>
                        end
                    end

                    % remove the objects
                    comps_keep = setdiff(this.mask_comps, comps_rm);
                    mask_obj = this.extract(comps_keep);
                otherwise
            end
        end

        function mask_obj = extract(this, r_)
            arguments
                this
                r_  (1,:)   % components or string label
            end

            if ~isstring(r_) && ~isnumeric(r_)
                throw(MException("regmask:invalidROI", ...
                    "ROI must be string or numeric vector."));
            end

            if isstring(r_)
                % transform string to label indices
                [~, loc] = ismember(r_, this.mask_lbl);
                indices = this.mask_comps(loc);
            else
                indices = r_;
                [~, loc] = ismember(r_, this.mask_comps);
            end

            % create new object
            lbl = this.MaskLabel(loc);
            mask = zeros(size(this.mask_vol), "like", this.mask_vol);
            for k = indices
                % extract components (indices)
                mask(this.mask_vol==k) = k;
            end

            mask_obj = regmask(mask, lbl);

        end

        function bdbox = getBoundingBox(this, r_)
            arguments
                this
                r_  (1,:)   % index or string label
            end

            if isstring(r_)
                [~, loc] = ismember(r_, this.mask_lbl);
            else
                [~, loc] = ismember(r_, this.mask_comps);
            end

            bdbox = this.bdbox(loc, :);
        end

        function save(this, folder, file)
            arguments
                this
                folder  (1,1)   string  {mustBeFolder}
                file    (1,1)   string
            end

            mask = this.MaskVol;
            rois = this.MaskROI;
            label = this.MaskLabel;

            % v7.3 need python hd5 support, save as v7 is convient
            save(fullfile(folder, file), "mask", "rois", "label", "-v7");
        end
    end

    methods(Access=private, Hidden)
        function calc_bounding_box(this)
            mask = this.mask_vol;
            comps = this.mask_comps;

            this.bdbox = nan(numel(comps), 6);

            for n = 1:numel(comps)
                bw_mask = (mask==comps(n));
                stats = regionprops3(bw_mask, "BoundingBox", "Volume");
                if isempty(stats), continue; end

                if size(stats, 1) > 1
                    [~, midx] = max(stats.Volume);
                    stats = stats(midx, :);
                end

                % generate kernel which shape as BoundingBox
                this.bdbox(n, :) = stats.BoundingBox; % (y,x,z)
            end
        end
    end

    methods(Static)
        function [mask_roi, mask_lbl] = vol2roi(mask_vol)
            % This function transforms volume representation to roi
            % representation
            arguments
                mask_vol    (:,:,:)
            end

            mask_roi = cell(size(mask_vol,3)+1, 1);
            mask_roi{1} = [size(mask_vol,[2,1,3]), 0];  % record the volume size
            lbl = [];

            for z = 1:size(mask_vol, 3)
                % region analysis on each plane
                mask_z = mask_vol(:,:,z);
                stats = regionprops("table", mask_z, "Area", "Centroid", ...
                    "EquivDiameter");
                if isempty(stats)
                    mask_roi{z+1} = double.empty(0, 4);
                    continue; 
                end

                lbl_z = find(stats.Area > 0);

                if ~isempty(lbl_z)
                    lbl = [lbl; lbl_z]; %#ok<AGROW>

                    mask_roi{z+1} = [stats.Centroid(lbl_z,:), stats.EquivDiameter(lbl_z,:)./2, ...
                        lbl_z];
                else
                    mask_roi{z+1} = double.empty(0, 4);
                end
            end

            if nargout == 2
                mask_lbl = reshape(string(unique(lbl, "sorted")), 1, []);
            end
        end

        function [mask_vol, mask_lbl] = roi2vol(mask_roi)
            % This function transforms roi representation to volume
            % representation
            arguments
                mask_roi    (:,1)   cell
            end

            sz = mask_roi{1}([2,1,3]);
            lbl = cellfun(@(x)x(:,4), mask_roi, "UniformOutput",false);
            label_max = max(cell2mat(lbl), [], "all");

            bits_list = 8*[1, 2, 4, 8];
            clist = nextpow2(label_max+1)./bits_list;
            dtype_bits = bits_list(find(clist<=1, 1, "first"));

            dtype = sprintf("uint%d", dtype_bits);
            mask_vol = zeros(sz, dtype);
            theta = reshape(linspace(0, 2*pi, 100), [], 1);
            unit_c = [cos(theta), sin(theta)];
            lbl = [];

            for z = 1:sz(3)
                roi_z = mask_roi{z+1};
                mask_z = zeros(sz(1:2), dtype);

                for k = 1:numel(roi_z)
                    % generate a circle and move the center
                    A = roi_z(k, 3)*unit_c + roi_z(k,1:2);

                    label_mask = roi_z(k,4)*poly2mask(A(:,1), A(:,2), sz(2), sz(1));
                    lbl = [lbl; roi_z(k,4)]; %#ok<AGROW>

                    mask_z = mask_z + cast(label_mask, dtype);
                end
                
                mask_vol(:,:,z) = mask_z;
            end

            if nargout == 2
                mask_lbl = reshape(string(unique(lbl, "sorted")), 1, []);
            end

        end

        function [mask_flat, mask_lbl] = flatten(maskobj)
            % This function transforms dense mask to sparse representation,
            % which could be recognize by python:scipy
            arguments
                maskobj (1,1)   regmask
            end

            mask = maskobj.MaskVol;

            d = unique(mask);
            d(d==0) = [];
            d(isnan(d)) = [];

            for k = 1:numel(d)
                vd_loc = find(mask == d(k));   % F order, linear sub index
                row_mask = [row_mask; vd_loc]; %#ok<AGROW>
                col_mask = [col_mask; k*ones(size(vd_loc))]; %#ok<AGROW>
            end
            vol_mask = ones(size(row_mask));

            % generate sparse matrix
            mask_flat = sparse(row_mask, col_mask, vol_mask, numel(mask), numel(d));

            mask_lbl = maskobj.MaskLabel;
        end

        function maskobj = recover(mask_flat, sz, mask_lbl)
            % This function transforms sparse mask to dense representation
            % where mask could be recognize by python:scipy
            arguments
                mask_flat   (:,:)   double  % sparse double
                sz          (1,3)   double  {mustBePositive, mustBeInteger}
                mask_lbl    (1,:)   string = ""
            end

            if ~issparse(mask_flat)
                throw(MException("regmask:recover_dense:invalidMaskFormat", ...
                    "Input mask must be sparse matrix."));
            end
            if size(mask_flat, 1) ~= prod(sz)
                throw(MException("regmask:recover_dense:invalidMaskSize", ...
                    "Invalid mask size."))
            end

            % reconstruct the mask volume
            v_mask = zeros(sz);
            for k = 1:size(mask_flat, 2)
                vidx = full(mask_flat(:, k));
                v_mask(vidx>0) = k;
            end

            maskobj = regmask(v_mask, mask_lbl);
        end

        function tf = overlap_gt(range_xy, bdbox, th)
            % This function calculates the intercept ratio for each bounding
            % box, if ratio is greater than th, return true
            arguments
                range_xy    (2,2)   double
                bdbox       (:,6)   double
                th          (:,1)   double  = 0.5
            end
            A = range_xy*[1, -1; 0, 1];
            A(:,2) = A(:,2) + 1;
            A = reshape(A, 1, []);
            B = round(bdbox(:,[1,2,4,5]));
            area = rectint(B, A);
            
            tf = (area./prod(bdbox(:,[4,5]), 2) > th);
        end
    end
end

function mustBeMask(A)
% A could be label matrix or sliced-circles cell
if ~isnumeric(A) && ~iscell(A)
    throw(MException("mustBeMask:invalidInput", ...
        "Input must be a label matrix or a sliced-circles cell."));
end

if isnumeric(A)
    if ~issparse(A) && (ndims(A)~=3 || any(A < 0, "all") || any(round(A)~=A, "all"))
        throw(MException("mustBeMask:invalidInput", ...
            "Input matrix must be a volume with nonnegtive integer or sparse matrix."));
    end
else
    % validate items
    if ~isvector(A) || any(cellfun(@(x)(size(x,2)~=4 || ~isnumeric(x)), A, "UniformOutput",true),"all")
        throw(MException("mustBeMask:invalidInput", ...
            "Input sliced-circles cell must be a vector with n-by-4 numeric matrix as items."));
    end
    if any(cellfun(@(x)(~isempty(x) && (any(round(x(:,4))~=x(:,4)) || any(x(:,4) < 0))), ...
            A, "UniformOutput",true))
        throw(MException("mustBeMask:invalidInput", ...
            "Circle indicator index must be nonnegtive integer."));
    end
end

end