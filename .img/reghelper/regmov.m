classdef regmov < matlab.mixin.Copyable
    %REGMOV This class is regmov definition, which construct the key data
    %struct for Register

    properties(Constant, Hidden)
        BYTES_UINT8 = 1
        BYTES_UINT16 = 2
        BYTES_UINT32 = 4
        BYTES_UINT64 = 8
        INNER_BLOCK_SIZE = 50
    end

    properties(Access=private, Hidden)
        mopt            % 1-by-12 table, with movie options
        t               % t-by-1 double, camera time
        zprojalg        % 1-by-1 string, method of z-projection
        zprojupd        % 1-by-1 logical, the flag for projection updating
        tform           % t-by-3 cell, with {global, local, manual} transformation unit
    end

    properties(Access=private, Hidden, NonCopyable)
        mptr            % 1-by-1 mpimg/mpimgs object or numeric array
        zdptr           % m-by-n-by(-c)-by-t uint16 array, z-projected data
    end

    properties(Access=public, Dependent)
        Bytes           % variable, get/___, not stored
        EMax            % variable, get/___, not stored
        EMedian         % variable, get/___, not stored
        EMin            % variable, get/___, not stored
        ETmplIdx        % variable, get/___, not stored
        Location        % variable, get/___, not stored
        MC              % variable, get/___, not stored
        MetaData        % variable, get/set, stored
        Movie           % variable, get/set, stored
        MovieZProj      % variable, get/___, stored
        RetainCache     % variable, get/set, not stored
        Time            % variable, get/set, stored
        Transformation  % variable, get/set, stored
        ZProjMethod     % variable, get/set, stored
    end

    methods
        function this = regmov(mptr_, mopt_, t_)
            %REGMOV A constructor
            arguments
                mptr_   
                mopt_   (1,12)  table  {mustBeMovieOptions}
                t_      (:,1)   double {mustBeNonnegative}
            end
            if ~ismember(class(mptr_), ["mpimg", "mpimgs"]) ...
                    && ~isnumeric(mptr_)
                throw(MException("regmov:invalidMoviePointer", ...
                    "Invalid movie pointer, Only mpimg, mpimgs object or array are supported."));
            else
                if ismember(class(mptr_), ["mpimg", "mpimgs"]) ...
                        && (numel(mptr_)  ~= 1)
                    throw(MException("regmov:invalidMoviePointer", ...
                        "Only one movie pointer is supported."));
                end
            end
            if numel(t_) ~= mopt_.frames
                throw(MException("regmov:badTime", ...
                    "Number of Timesteps not match."));
            end

            this.mptr = mptr_;
            this.mopt = mopt_;
            this.t = t_;
            this.zprojalg = "max";
            this.tform = cell(mopt_.frames, 3);
            this.zprojupd = false;
        end

        function r = get.Movie(this)
            % this operation will call deep copy (constructor) if mptr is
            % mpimg or mpimgs object
            if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                % hide other properties
                r = this.mptr.Data;
            elseif isnumeric(this.mptr)
                r = this.mptr;
            end
        end

        function set.Movie(this, r_)
            arguments
                this
                r_  
            end

            if ~ismember(class(r_), ["mpimg", "mpimgs"]) ...
                    && ~isnumeric(r_)
                throw(MException("regmov:invalidMoviePointer", ...
                    "Invalid movie pointer, Only mpimg, mpimgs object or array are supported."));
            else
                if ismember(class(r_), ["mpimg", "mpimgs"]) ...
                        && (numel(r_)  ~= 1)
                    throw(MException("regmov:invalidMoviePointer", ...
                        "Only one movie pointer is supported."));
                end
            end

            if (ismember(class(this.mptr), ["mpimg", "mpimgs"]) ...
                    && ismember(class(r_), ["mpimg", "mpimgs"])) ...
                    || (isnumeric(this.mptr) && isnumeric(r_))
                this.mptr = r_;
            elseif ismember(class(this.mptr), ["mpimg", "mpimgs"]) ...
                    && isnumeric(r_)
                this.mptr.Data = r_;
            elseif isnumeric(this.mptr) ...
                    && ismember(class(r_), ["mpimg", "mpimgs"])
                this.mptr = r_.Data;
            end
        end

        function r = get.MetaData(this)
            r = this.mopt;
        end

        function set.MetaData(this, r_)
            arguments
                this
                r_  (1,12)  table   {mustBeMovieOptions}
            end

            this.mopt = r_;
        end

        function r = get.Time(this)
            r = this.t;
        end

        function set.Time(this, r_)
            arguments
                this
                r_  (:,1)   double {mustBeNonnegative}
            end

            if numel(r_) ~= this.mopt.frames
                throw(MException("regmov:badTime", ...
                    "Number of Timesteps not match."));
            end

            this.t = r_;
        end

        function r = get.MovieZProj(this)
            % NOTE: This function is not real-time function, do not call in
            % real-time situation

            if isempty(this.zdptr)
                % update z projection source
                update_flag = true;
            else
                update_flag = this.zprojupd;
            end

            if update_flag == true
                if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                    [vz, pz] = ismember("Z", this.mptr.DimOrder);
                    if ~vz
                        throw(MException("regmov:invalidOperation", ...
                            "Movie does not have dimension: 'Z'."));
                    end

                    rsp_t = repmat({':'}, 1, this.mptr.DataDims);
                    rsp_tnz = repmat({':'}, 1, this.mptr.DataDims-1);
                    [~, tloc] = ismember("T", this.mopt.dimOrder);
                    [~, zloc] = ismember("Z", this.mopt.dimOrder);

                    % partial loading and calculating
                    block_n = ceil(this.mopt.frames/this.INNER_BLOCK_SIZE);
                    mov_sz_zproj = this.mptr.DataSize;
                    mov_sz_zproj(zloc) = [];
                    r = zeros(mov_sz_zproj, this.mptr.DataType);

                    % slow 50% than exclude typically
                    for nb = 1:block_n
                        vidx = ((nb-1)*this.INNER_BLOCK_SIZE+1):min(...
                            nb*this.INNER_BLOCK_SIZE, this.mopt.frames);
                        rsp_t{tloc} = vidx;
                        rsp_tnz{tloc-(tloc>zloc)} = vidx;

                        data = subsref(this.mptr.Data, struct('type','()','subs',{rsp_t}));

                        dr = Projection(data, this.zprojalg, pz);

                        r = subsasgn(r, struct('type','()','subs',{rsp_tnz}), dr);
                    end
                elseif isnumeric(this.mptr)
                    [vz, pz] = ismember("Z", this.mopt.dimOrder);
                    if ~vz
                        throw(MException("regmov:invalidOperation", ...
                            "Movie does not have dimension: 'Z'."));
                    end

                    r = Projection(this.mptr, this.zprojalg, pz);
                end

                % update
                this.zdptr = r;
            else
                % copy
                r = this.zdptr;
            end
        end

        function r = get.ZProjMethod(this)
            r = this.zprojalg;
        end

        function set.ZProjMethod(this, r_)
            arguments
                this
                r_  (1,1)   string {mustBeMember(r_, ["max","min", ...
                    "median","mean"])} = "max"
            end

            if ~isequal(r_, this.zprojalg)
                % update projection data
                this.zprojalg = r_;

                this.zprojupd = true;   % force to update
                this.MovieZProj;    % 
                this.zprojupd = false;  % reset
            end
        end

        function r = get.Transformation(this)
            r = this.tform;
        end

        function set.Transformation(this, r_)
            arguments
                this
                r_  (:, 3) cell
            end
            if size(r_, 1) ~= this.mopt.frames
                throw(MException("regmov:badTransformation", ...
                    "Number of Timesteps not match."));
            end

            this.tform = r_;
        end

        function r = get.Location(this)
            if ~ismember(class(this.mptr), ["mpimg", "mpimgs"])
                r = "Memory";
            else
                r = string(this.mptr.FileName);
            end
        end

        function r = get.Bytes(this)
            mopt_ = this.mopt; %#ok<NASGU>
            t_ = this.t; %#ok<NASGU>
            zdptr_ = this.zdptr; %#ok<NASGU>
            tform_ = this.tform; %#ok<NASGU>
            prop_ = struct2table(whos("mopt_","t_","zdptr_","tform_"));
            r.mem = sum(prop_.bytes);

            if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                r.map = this.mptr.DataBytes;
            elseif isnumeric(this.mptr)
                r.map = 0;
                switch class(this.mptr)
                    case {'uint8', 'int8'}
                        bytes_per_elem = this.BYTES_UINT8;
                    case {'uint16','int16'}
                        bytes_per_elem = this.BYTES_UINT16;
                    case {'uint32','int32','single'}
                        bytes_per_elem = this.BYTES_UINT32;
                    case {'uint64','int64','double'}
                        bytes_per_elem = this.BYTES_UINT64;
                end
                r.mem = r.mem + numel(this.mptr)*bytes_per_elem;
            end
        end

        function r = get.EMin(this)
            % use random re-sample for quick and better lower bound 
            % estimation

            r = zeros([1, this.MetaData.channels]);

            npix = this.MetaData.width * this.MetaData.height * ...
                this.MetaData.images;

            % we select at most 1e6 points from
            LX = rand(min(1e6, npix), 5);
            LX(:,1) = round((this.MetaData.height-1)*LX(:,1)) + 1;  % rows
            LX(:,2) = round((this.MetaData.width-1)*LX(:,2))+ 1;   % cols
            LX(:,4) = round((this.MetaData.slices-1)*LX(:,4)) + 1;  % depth
            LX(:,5) = round((this.MetaData.frames-1)*LX(:,5)) + 1;  % time

            % modify the color channel indices
            for k = 1:numel(r)
                LX(:,3) = k;
                ind_k = sub2ind(size(this.Movie), LX(:,1),LX(:,2),LX(:,3),LX(:,4),LX(:,5));
                r(k) = min(this.Movie(ind_k),[],"all");
            end
        end

        function r = get.EMax(this)
            % use frame sorted re-sample for quick and better upper
            % bound estimation

            r = zeros([1, this.MetaData.channels]);

            floc = randperm(this.MetaData.frames, 5);

            % modify the color channel indices
            for k = 1:numel(r)
                r(k) = max(this.Movie(:,:,k,:,floc), [], "all");    % dimensionalities unsafe
            end
        end

        function r = get.EMedian(this)
            % use frame sorted re-sample for quick and better median
            % estimation

            r = zeros([1, this.MetaData.channels]);

            floc = randperm(this.MetaData.frames, 5);

            % modify the color channel indices
            for k = 1:numel(r)
                r(k) = median(this.Movie(:,:,k,:,floc),"all");      % dimensionalities unsafe
            end
        end

        function r = get.ETmplIdx(this)
            % use summed gradient on median plane
            % note that aligned volume must be warped, interpolation 
            % caused blur so that the gradient sum should be smaller 
            % than template volume

            % sobel operator
            sobel_x = [-1 0 1; -2 0 2; -1 0 1];
            sobel_y = [-1 -2 -1; 0 0 0; 1 2 1];

            % only use the red channel
            cidx = (this.mopt.cOrder == "r");
            if any(cidx)
                sidies = ceil(linspace(1, this.mopt.slices, min(this.mopt.slices, 7)));
                % throw first and last slice, keep 5 slices at most
                sidies = sidies(2:end-1);
                sumGrads = zeros(this.mopt.frames, 1);

                % for loop is faster than 4-D operation because memory
                % reallocating, speed up about 40%
                for sidx = sidies
                    % resample some planes
                    imgs = im2double(squeeze(this.Movie(:,:,cidx,sidx,:)));     % dimensionalities unsafe

                    % calculate the sum of intensity gradient
                    Gx = imfilter(imgs, sobel_x, "replicate");
                    Gy = imfilter(imgs, sobel_y, "replicate");
                    sumGrads = sumGrads + squeeze(sum(sqrt(Gx.^2 + Gy.^2), [1,2]));  % keep time dimensionalities
                end

                % return the maxima location of summed gradient
                [~, r] = max(sumGrads);

                % maxima should also be an outlier
                if ~ismember(r, find(isoutlier(sumGrads)))
                    disp("The estimated template index may be not good.");
                end
            else
                r = 0;
            end
        end

        function r = get.MC(this)
            % this function get the mass center of each volume, 
            % r as t-by-3-by-c array, with [x,y,z]

            [~, ctloc] = ismember(["C","T"], this.mopt.dimOrder);
            dimsum = setdiff(1:numel(this.mopt.dimOrder), ctloc);
            rsp = repmat({1}, 1, numel(this.mopt.dimOrder));

            if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                % note that mpimg as low level abstract of image, no reason
                % to support mass-center
                % we calculate mass center of volumes and concatenate them
                xloc = find("X"==this.mopt.dimOrder); ny = this.mptr.DataSize(xloc);
                yloc = find("Y"==this.mopt.dimOrder); nx = this.mptr.DataSize(yloc);
                zloc = find("Z"==this.mopt.dimOrder); nz = this.mptr.DataSize(zloc);
                rsp_x = rsp; rsp_x{yloc} = [];  xgrid = reshape(1:nx, rsp_x{:});
                rsp_y = rsp; rsp_y{xloc} = [];  ygrid = reshape(1:ny, rsp_y{:});
                rsp_z = rsp; rsp_z{zloc} = [];  zgrid = reshape(1:nz, rsp_z{:});
                mc_ref_x = zeros(2, this.mopt.frames);
                mc_ref_y = zeros(2, this.mopt.frames);
                mc_ref_z = zeros(2, this.mopt.frames);
                
                cloc = ctloc(1); tloc = ctloc(2);
                rsp_t = repmat({':'}, 1, this.mptr.DataDims);
                rsp_c = rsp; rsp_c{cloc} = []; 
                cgrid = uint16(reshape(numel(this.mopt.cOrder)*constdef.CAMERA_BACKGROUND, rsp_c{:}));

                % partial loading and calculating
                block_n = ceil(this.mopt.frames/this.INNER_BLOCK_SIZE);
                
                % slow 40% than exclude typically
                for nb = 1:block_n
                    vidx = ((nb-1)*this.INNER_BLOCK_SIZE+1):min(...
                        nb*this.INNER_BLOCK_SIZE, this.mopt.frames);
                    rsp_t{tloc} = vidx;

                    data = subsref(this.mptr.Data, struct('type','()','subs',{rsp_t}));
                    pixval_tot = double(squeeze(sum(data-cgrid, dimsum)));

                    % calculate mass center
                    mc_ref_x(:, vidx) = reshape(sum(xgrid.*sum(data-cgrid, [xloc,zloc]), yloc), ...
                        [this.mptr.DataSize(ctloc(1)), numel(vidx)])./(pixval_tot+eps);
                    mc_ref_y(:, vidx) = reshape(sum(ygrid.*sum(data-cgrid, [yloc,zloc]), xloc), ...
                        [this.mptr.DataSize(ctloc(1)), numel(vidx)])./(pixval_tot+eps);
                    mc_ref_z(:, vidx) = reshape(sum(zgrid.*sum(data-cgrid, [xloc,yloc]), zloc), ...
                        [this.mptr.DataSize(ctloc(1)), numel(vidx)])./(pixval_tot+eps);
                end
            elseif isnumeric(this.mptr)
                % get X,Y,Z dimension location and generate grid
                xloc = find("X"==this.mopt.dimOrder); ny = size(this.mptr, xloc);
                yloc = find("Y"==this.mopt.dimOrder); nx = size(this.mptr, yloc);
                zloc = find("Z"==this.mopt.dimOrder); nz = size(this.mptr, zloc);
                rsp_x = rsp; rsp_x{yloc} = [];  xgrid = reshape(1:nx, rsp_x{:});
                rsp_y = rsp; rsp_y{xloc} = [];  ygrid = reshape(1:ny, rsp_y{:});
                rsp_z = rsp; rsp_z{zloc} = [];  zgrid = reshape(1:nz, rsp_z{:});
                rsp_c = rsp; rsp_c{ctloc(1)} = []; 
                cgrid = uint16(reshape(numel(this.mopt.cOrder)*constdef.CAMERA_BACKGROUND, rsp_c{:}));
                pixval_tot = double(squeeze(sum(this.mptr-cgrid, dimsum)));
                pixval_tot = reshape(pixval_tot, numel(cgrid), []);

                % calculate mass center
                mc_ref_x = reshape(sum(xgrid.*sum(this.mptr-cgrid, [xloc,zloc]), yloc), ...
                    size(this.mptr, ctloc))./(pixval_tot+eps);
                mc_ref_y = reshape(sum(ygrid.*sum(this.mptr-cgrid, [yloc,zloc]), xloc), ...
                    size(this.mptr, ctloc))./(pixval_tot+eps);
                mc_ref_z = reshape(sum(zgrid.*sum(this.mptr-cgrid, [xloc,yloc]), zloc), ...
                 size(this.mptr, ctloc))./(pixval_tot+eps);
            end

            r = permute(cat(3, mc_ref_x, mc_ref_y, mc_ref_z), [2,3,1]);
        end

        function r = get.RetainCache(this)
            if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                r = this.mptr.RetainCache;
            else
                r = false;  % auto clean by MATLAB garbage collector
            end
        end

        function set.RetainCache(this, r)
            arguments
                this 
                r       (1,1)   logical
            end

            if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                this.mptr.RetainCache = r;
            else
                % do nothing
            end
        end
    end

    methods(Access = protected)
        % Override copyElement method:
        function cpt = copyElement(this)
            % Make a shallow copy of all properties except memptr
            cpt = copyElement@matlab.mixin.Copyable(this);

            if isnumeric(this.mptr)
                % allocate after changing, MATLAB backend takes over,
                % implicit controller
                cpt.mptr = this.mptr;   
            elseif ismember(class(this.mptr), ["mpimg", "mpimgs"])
                % allocate right now, programming takes over manually
                % explicit controller
                cpt.mptr = this.mptr.copy();
            end
        end
    end

    methods(Access=public, Hidden)
        function gather(this)
            % this function gathers data from disk file to memory
            if isnumeric(this.mptr)
                % already in memory
                return;
            else
                D = this.mptr.Data;
                delete(this.mptr);  % release resource
                this.mptr = D;
            end
        end

        function spread(this)
            % this function spreads data from memory to disk file
            if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                % already on disk
                return;
            else
                tmpfolder = mpimg.findtmpfolder(this.mopt); % mpimgs ?
                D = this.mptr;
                this.mptr = [];     % free variable
                this.mptr = mpimg(tmpfolder, [], D, this.mopt.dimOrder);
            end
        end
    end

    methods(Access=public)
        function movobj = vcrop(this, dim_, r_)
            arguments
                this
                dim_  (1,1) string
                r_    (:,2) double {mustBePositive, mustBeInteger}
            end

            validateattributes(r_', "double", "nondecreasing");

            mopt_ = this.mopt;
            t_ = this.t;
            tf_ = this.tform;

            if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                switch dim_
                    case "XY"
                        mptr_ = this.mptr.vcrop("YX", r_);
                        mopt_.width = diff(r_(1,:))+1;
                        mopt_.height = diff(r_(2,:))+1;
                        % crop the non-rigid registration displacement field
                        tf_(:,2) = cellfun(@(x)df_crop_xy(x, r_), ...
                            tf_(:,2),'UniformOutput',false);
                    case "Z"
                        mptr_ = this.mptr.vcrop("Z", r_);
                        mopt_.slices = diff(r_)+1;
                        % crop the non-rigid registration displacement field
                        tf_(:,2) = cellfun(@(x)df_crop_z(x, r_), ...
                            tf_(:,2),'UniformOutput',false);
                    otherwise
                        throw(MException("regmov:vcrop:invalidUse", ...
                            "Undefined using of vcrop. Only 'XY' and 'Z' mode are supported."));
                end

                % update images number (if slices changed)
                mopt_.images = mopt_.slices*mopt_.channels*mopt_.frames;

                % generate a regmov obj
                movobj = regmov(mptr_, mopt_, t_);
                movobj.Transformation = tf_;
            elseif isnumeric(this.mptr)
                switch dim_
                    case "XY"
                        mov = CropXY(this.Movie, r_');
                        mopt_.width = diff(r_(1,:))+1;
                        mopt_.height = diff(r_(2,:))+1;
                        tf_(:,2) = cellfun(@(x)df_crop_xy(x, r_), tf_(:,2),'UniformOutput',false);
                    case "Z"
                        mov = CropZ(this.Movie, r_);
                        mopt_.slices = diff(r_)+1;
                        tf_(:,2) = cellfun(@(x)df_crop_z(x, r_), tf_(:,2),'UniformOutput',false);
                    otherwise
                        throw(MException("regmov:vcrop:invalidUse", ...
                            "Undefined using of vcrop. Only 'XY' and 'Z' mode are supported."));
                end
                mopt_.images = mopt_.slices*mopt_.channels*mopt_.frames;

                movobj = regmov(mov, mopt_, t_);
                movobj.Transformation = tf_;
            end

            function A = df_crop_xy(A, r)
                % A: m-by-n-by-p-by-3
                if ~isempty(A)
                    A = A(r(2,1):r(2,2),r(1,1):r(1,2),:,:);
                end
            end

            function A = df_crop_z(A, r)
                % A: m-by-n-by-p-by-3
                if ~isempty(A)
                    A = A(:,:,r(1):r(2),:);
                end
            end
        end

        function movobj = tcrop(this, r_)
            arguments
                this
                r_    (1,:) double {mustBePositive, mustBeInteger}
            end

            validateattributes(r_, "double", "nondecreasing");

            mopt_ = this.mopt; 

            mopt_.frames = numel(r_);
            mopt_.images = mopt_.slices*mopt_.channels*mopt_.frames;
            t_ = this.t(r_);
            tf_ = this.tform(r_, :);

            if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                mptr_ = this.mptr.tcrop(r_);
                movobj = regmov(mptr_, mopt_, t_);
            else
                mov = CropT(this.mptr, r_);
                movobj = regmov(mov, mopt_, t_);
            end

            movobj.Transformation = tf_;
        end

        function ctset(this, v, cr, tr)
            % This function set v to raw data with location: c & t dimension
            arguments
                this
                v   (:,:,:,:)
                cr  (1,1)    string {mustBeMember(cr, ["r","g","b"])}
                tr  (1,1)    string
            end

            cr = (cr == this.mopt.cOrder);
            tr = str2num(tr); %#ok<ST2NM>

            if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                this.mptr.ctset(v, cr, tr);
            else
                % calling inner array processing
                this.ctset_(v, cr, tr);
            end
        end

        function fliplr(this)
            % This function flips left and right on XY plane
            if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                this.mptr.fliplr();
            else
                % calling inner array processing
                this.mptr = fliplr(this.mptr);
            end
        end

        function flipud(this)
            % This function flips up and down on XY plane
            if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                this.mptr.flipud();
            else
                % calling inner array processing
                this.mptr = flipud(this.mptr);
            end
        end

        function rot90(this, k)
             % This function rotate 90 degrees on XY plane k times, 1 as default
             arguments
                 this
                 k  (1,1)   double {mustBePositive,mustBeInteger} = 1
             end

             if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                 this.mptr.rot90(k);
             else
                 % calling inner array processing
                 this.mptr = rot90(this.mptr, k);
             end
             
             % odd rotation lead to X&Y size exchanging
             if mod(k, 2)~=0
                 width_old = this.mopt.width;
                 xRes_old = this.mopt.xRes;
                 this.mopt.width = this.mopt.height;
                 this.mopt.xRes = this.mopt.yRes;
                 this.mopt.height = width_old;
                 this.mopt.yRes = xRes_old;
             end
        end

        function movobj = maskby(this, mask)
            % This function masked data and return a new regmov object
            arguments
                this    (1,1)
                mask    (1,1)   regmask
            end

            mkvol = mask.MaskVol;

            if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                this.mptr.maskby(mkvol);
            else
                this.mptr = MaskMovie(this.mptr, mkvol);
            end

            movobj = this;
        end

        function delete(this)
            if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                delete(this.mptr);
            else
                % free memory
                this.mptr = [];
            end

            % ~
            clear this
        end

        function r = isempty(this)
            r = isempty(this.Movie);
        end
    end

    methods(Access=private, Hidden)
        function expr = gen_expstr(this, cr, tr, istr)
            t_loc = find(this.MetaData.dimOrder=="T");
            c_loc = find(this.MetaData.dimOrder=="C");
            cr = string(find(cr));

            if isempty(cr) || isempty(t_loc) || isempty(c_loc)
                throw(MException("regmov:ctset:invalidMovie", ...
                    "Bad calling: invalid movie dimension."));
            end

            mov_ndim = numel(this.MetaData.dimOrder);

            expr = "";

            for dp = 1:mov_ndim
                if dp == c_loc
                    expr = expr + cr;
                elseif dp == t_loc
                    expr = expr + tr;
                else
                    expr = expr + ":";
                end
                if dp ~= mov_ndim, expr = expr + ","; end
            end

            expr = sprintf("this.mptr(%s)=%s;", expr, istr);
        end

        function ctset_(this, v, cr, tr)
            % optimize the usual case for fast saving
            switch this.mopt.dimOrder.join("")
                case "XYCZT"
                    this.mptr(:,:,cr,:,tr) = v;
                case "XYZCT"
                    this.mptr(:,:,:,cr,tr) = v;
                otherwise
                    % generate the crop expression
                    expr = this.gen_expstr(cr, tr, "v");

                    % simple modifying
                    eval(expr);
            end
        end
    end

    methods(Static)
        % This function create an empty regmov object
        function mov = empty()
            mopt_ = table('Size',[1,12], 'VariableTypes',{'double','double',...
                'double','double','double','double','double','double','double',...
                'string','string','string'},'VariableNames',{'width','height',...
                'slices','channels','frames','images','xRes','yRes','zRes',...
                'dataType','dimOrder','cOrder'});
            mopt_.dataType = "uint16";
            mopt_.dimOrder = string.empty(1,0);
            mopt_.cOrder = string.empty(1,0);

            mov = regmov(uint16.empty(), mopt_, []);
        end
    end
end


% ====================== local utility function ===================
function mustBeMovieOptions(A)
MOVIE_OPTIONS_FIELD = ["width", "height", "slices", "channels", "frames", ...
    "images", "xRes", "yRes", "zRes", "dataType","dimOrder", "cOrder"];
varnames = A.Properties.VariableNames;
if ~isempty(setxor(varnames, MOVIE_OPTIONS_FIELD))
    throwAsCaller(...
        createValidatorExceptionWithValue(createPrintableList(MOVIE_OPTIONS_FIELD), ...
        'regmov:validators:mustBeEqualGenericText',...
        'regmov:validators:mustBeMember')...
        );
end

validateattributes(A.width, "double", {'scalar','integer','nonnegative'});
validateattributes(A.height, "double", {'scalar','integer','nonnegative'});
validateattributes(A.slices, "double", {'scalar','integer','nonnegative'});
validateattributes(A.channels, "double", {'scalar','integer','nonnegative'});
validateattributes(A.frames, "double", {'scalar','integer','nonnegative'});
validateattributes(A.images, "double", {'scalar','integer','nonnegative'});
validateattributes(A.xRes, "double", {'scalar','nonnegative'});
validateattributes(A.yRes, "double", {'scalar','nonnegative'});
validateattributes(A.zRes, "double", {'scalar','nonnegative'});
validateattributes(A.dataType, "string", {'scalar'});
validateattributes(A.dimOrder, "string", {'vector'});
validateattributes(A.cOrder, "string", {'vector'});
end
