classdef regmask < handle
    %REGMASK This class is mask object defination, with basic operation on
    %mask
    
    properties(Access=private, Hidden)
        mask_vol
        mask_roi
        bdbox
        label
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

        end

        function mask_vol = roi2vol(mask_roi)

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