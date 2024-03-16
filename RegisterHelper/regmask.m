classdef regmask < handle
    %REGMASK This class is mask object defination, with basic operation on
    %mask
    
    properties(Access=private, Hidden)
        mask_vol        % m-by-n-by-p label volume
        mask_roi        % p+1-by-1 cell array, each item for one plane, first item with background
        bdbox           % 1-by-s cell, each item is components bounding box
        label           % 1-by-s string array, each item is a component label
    end

    properties(GetAccess=public, Dependent)
        MaskVol
        MaskROI
        MaskLabel
    end
    
    methods
        function this = regmask(mask_, label_)
            %REGMASK A constructor
            arguments
                mask_   {mustBeMask}
                label_  (1,:)   string
            end

            % transform mask and stored
            if iscell(mask_)
                this.mask_roi = mask_;
                this.mask_vol = regmask.roi2vol(mask_);
            else
                this.mask_vol = mask_;
                this.mask_roi = regmask.vol2roi(mask_);
            end

            % 
            this.label = label_;

            %
            this.calc_bounding_box();
        end

        function r = get.MaskVol(this)
            r = this.mask_vol;
        end
        
        function r = get.MaskROI(this)
            r = this.mask_roi;
        end

        function r = get.MaskLabel(this)
            r = this.label;
        end

        function crop(this, r_)
            arguments
                this
                r_   (:,2)   {mustBePositive, mustBeInteger}
            end
        end

        function maskobj = extractBefore(this, r_)
            arguments
                this
                r_  (1,1)   % index or label
            end


        end

        function maskobj = extractAfter(this, r_)
            arguments
                this
                r_  (1,1)   % index or label
            end

        end

        function maskobj = extractBetween(this, r1_, r2_)
            arguments
                this
                r1_  (1,1)  % index or label
                r2_  (1,1)  % index or label
            end

        end

        function bdbox = getBoundingBox(this, r_)
            arguments
                this
                r_  (1,1)   % index or label
            end


        end
    end

    methods(Access=private, Hidden)
        function calc_bounding_box(this)

        end
    end

    methods(Static, Hidden)
        function mask_roi = vol2roi(mask_vol)
            % This function transform volume representation to roi
            % representation
            arguments
                mask_vol    (:,:,:)
            end

            mask_roi = cell(size(mask_vol,3)+1, 1);
            mask_roi{1} = [size(mask_roi,[2,1,3]), 0];  % record the volume size

            for z = 1:size(mask_vol, 3)
                % region analysis on each plane
                mask_z = mask_vol(:,:,z);
                stats = regionprops("table", mask_z, "Volume", "Centroid", ...
                    "EquivDiameter","PixelIdxList");
                stats(stats.Area==0) = [];
                if ~isempty(stats)
                    pidx = cellfun(@(x)(x(round(numel(x)/2))), stats.PixelIdxList, ...
                        "UniformOutput",true);
                    lbl_z = mask_z(pidx);

                    mask_roi{z+1} = [stats.Centroid, stats.EquivDiameter./2, ...
                        lbl_z];
                else
                    mask_roi{z+1} = double.empty(0, 4);
                end
            end
        end

        function mask_vol = roi2vol(mask_roi)
            arguments
                mask_roi    (:,1)   cell
            end

            sz = mask_roi{1}([2,1,3]);
            label_max = max(cellfun(@(x)max(x(:,4)), mask_roi, ...
                "UniformOutput",true), [], "all");

            bits_list = 8*[1, 2, 4, 8];
            clist = nextpow2(label_max+1)./bits_list;
            dtype_bits = bits_list(find(clist<=1, 1, "first"));

            dtype = sprintf("uint%d", dtype_bits);
            mask_vol = zeros(sz, dtype);
            theta = reshape(linspace(0, 2*pi, 100), [], 1);
            unit_c = [cos(theta), sin(theta)];

            for z = 1:sz(3)
                roi_z = mask_roi{z+1};
                mask_z = zeros(sz(1:2), dtype);

                for k = 1:numel(roi_z)
                    % generate a circle and move the center
                    A = roi_z(k, 3)*unit_c + roi_z(k,1:2);

                    label_mask = roi_z(k,4)*poly2mask(A(:,1), A(:,2), sz(2), sz(1));

                    mask_z = mask_z + cast(label_mask, dtype);
                end
                
                mask_vol(:,:,z) = mask_z;
            end
        end

        function mask_flat = sparse_flatten(maskobj)

        end

        function maskobj = recover_dense(mask_flat, sz)

        end
    end
end

function mustBeMask(A)
% A could be label matrix or sliced-circles cell
if ~ismatrix(A) && ~iscell(A)
    throw(MException("mustBeMask:invalidInput", ...
        "Input must be a label matrix or a sliced-circles cell."));
end

if ismatrix(A)
    if ~issparse(A) && (ndims(A)~=3 || ~isnumeric(A) || any(A < 0) || ...
            any(round(A)~=A, "all"))
        throw(MException("mustBeMask:invalidInput", ...
            "Input matrix must be a volume with nonnegtive integer or sparse matrix."));
    end
else
    % validate items
    if ~isvector(A) || any(cellfun(@(x)(size(x,2)~=4 || ~isnumeric(x)), A, "UniformOutput",true))
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