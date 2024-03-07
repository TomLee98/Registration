classdef regmov < matlab.mixin.Copyable
    %REGMOV This class is regmov defination, which construct the key data
    %struct for Register

    properties(Constant, Hidden)
        BYTES_UINT8 = 1
        BYTES_UINT16 = 2
        BYTES_UINT32 = 4
        BYTES_UINT64 = 8
    end

    properties(Access=private, Hidden)
        mopt            % 1-by-12 table, with movie options
        t               % t-by-1 double, camera time
        zprojalg        % 1-by-1 string, method of z-projection
        tform           % t-by-3 cell, with {global, local, manual} transformation unit
        bkg             % c-by-1 double, c for color channel number
    end

    properties(Access=private, Hidden, NonCopyable)
        mptr            % 1-by-1 mpimg/mpimgs object or numeric array
    end

    properties(Access=public, Dependent)
        Movie           % variable, get/set, stored
        MetaData        % variable, get/set, stored
        Time            % variable, get/set, stored
        MovieZProj      % variable,     get, not stored
        ZProjMethod     % variable, get/set, stored
        Transformation  % variable, get/set, stored
        Bytes           % variable,     get, not stored
        Background      % variable,     get, stored
        EMin            % variable,     get, not stored
        EMax            % variable,     get, not stored
        EMedian         % variable,     get, not stored
    end

    methods
        function this = regmov(mptr_, mopt_, t_, bkg_)
            %REGMOV A constructor
            arguments
                mptr_   
                mopt_   (1,12)  table  {mustBeMovieOptions}
                t_      (:,1)   double {mustBeNonnegative}
                bkg_            double {mustBeNonnegative}
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
            if numel(bkg_) ~= numel(mopt_.cOrder)
                throw(MException("regmov:badTime", ...
                    "Number of color channels not match."));
            end

            this.mptr = mptr_;
            this.mopt = mopt_;
            this.t = t_;
            this.zprojalg = "max";
            this.tform = cell(mopt_.frames, 3);
            this.bkg = bkg_;
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

            if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                [vz, pz] = ismember("Z", this.mptr.DimOrder);
                if ~vz
                    throw(MException("regmov:invalidOperation", ...
                        "Movie does not have dimension: 'Z'."));
                end

                r = Projection(this.mptr.Data, this.zprojalg, pz);
            elseif isnumeric(this.mptr)
                [vz, pz] = ismember("Z", this.mopt.dimOrder);
                if ~vz
                    throw(MException("regmov:invalidOperation", ...
                        "Movie does not have dimension: 'Z'."));
                end

                r = Projection(this.mptr, this.zprojalg, pz);
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

            this.zprojalg = r_;
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

        function r = get.Bytes(this)
            mopt_ = this.mopt; %#ok<NASGU>
            t_ = this.t; %#ok<NASGU>
            zprojalg_ = this.zprojalg; %#ok<NASGU>
            tform_ = this.tform; %#ok<NASGU>
            bkg_ = this.bkg; %#ok<NASGU>
            prop_ = struct2table(whos("mopt_","t_","zprojalg_","tform_","bkg_"));
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

        function r = get.Background(this)
            r = this.bkg;
        end

        function r = get.EMin(this)
            % use random resample for quick and better lower bound 
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
            % use frame sortted resample for quick and better upper
            % bound estimation

            r = zeros([1, this.MetaData.channels]);

            floc = randperm(this.MetaData.frames, 5);

            % modify the color channel indices
            for k = 1:numel(r)
                r(k) = max(this.Movie(:,:,k,:,floc), [], "all");
            end
        end

        function r = get.EMedian(this)
            % use frame sortted resample for quick and better median
            % estimation

            r = zeros([1, this.MetaData.channels]);

            floc = randperm(this.MetaData.frames, 5);

            % modify the color channel indices
            for k = 1:numel(r)
                r(k) = median(this.Movie(:,:,k,:,floc),"all");
            end
        end
    end

    methods(Access = protected)
        % Override copyElement method:
        function cpt = copyElement(this)
            % Make a shallow copy of all properties except memptr
            cpt = copyElement@matlab.mixin.Copyable(this);

            if isnumeric(this.mptr)
                cpt.mptr = this.mptr;
            elseif ismember(class(this.mptr), ["mpimg", "mpimgs"])
                % Make a deep copy of the mpimg/mpimgs object
                cpt.mptr = this.mptr.copy();
            end
        end
    end

    methods(Access=public)
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
                % already in memory
                return;
            else
                tmpfolder = mpimg.findtmpfolder(this.mopt); % mpimgs ?
                D = this.mptr;
                this.mptr = [];     % free variable
                this.mptr = mpimg(tmpfolder, [], D, this.mopt.dimOrder);
            end
        end

        function delete(this)
            if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                delete(this.mptr);
            end

            % ~
        end
    end

    methods(Access=public, Hidden)
        function r = isempty(this)
            r = isempty(this.Movie);
        end

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
            bkg_ = this.bkg;

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
                movobj = regmov(mptr_, mopt_, t_, bkg_);
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

                movobj = regmov(mov, mopt_, t_, this.bkg);
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
                movobj = regmov(mptr_, mopt_, t_, this.bkg);
            else
                mov = CropT(this.mptr, r_);
                movobj = regmov(mov, mopt_, t_, this.bkg);
            end

            movobj.Transformation = tf_;
        end

        function rmbkg(this, r_)
            % This function remove the background and store it
            arguments
                this
                r_  double  {mustBeNonnegative, mustBeVector}
            end

            if numel(r_) ~= this.MetaData.channels
                throw(MException("regmov:rmbkg:invalidChannelsNumber", ...
                    "Channels number not match."));
            end

            this.bkg = r_;

            if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                % mpimg overloading 'minus' operator
                this.mptr = this.mptr - r_;
            else
                % remove on each channel
                locstr = repmat(":", 1, 5);
                for k = 1:numel(r_)
                    locstr("C"==this.MetaData.dimOrder) = string(k);
                    d_expstr = "this.mptr(" + locstr.join(",") + ")";
                    sub_str = sprintf("-this.bkg(%d);", k);
                    eval(d_expstr + "=" + d_expstr + sub_str);
                end
            end
        end

        function rcbkg(this)
            % This function recover the background and reset it to zero
            if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                this.mptr = this.mptr + this.bkg;
            else
                 % add on each channel
                locstr = repmat(":", 1, 5);
                for k = 1:numel(this.bkg)
                    locstr("C"==this.MetaData.dimOrder) = string(k);
                    d_expstr = "this.mptr(" + locstr.join(",") + ")";
                    sub_str = sprintf("+this.bkg(%d);", k);
                    eval(d_expstr + "=" + d_expstr + sub_str);
                end
            end

            this.bkg = 0*this.bkg;
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

                    % simple modifing
                    eval(expr);
            end
        end
    end

    methods(Static, Hidden)
        function mov = empty()
            mopt_ = table('Size',[1,12], 'VariableTypes',{'double','double',...
                'double','double','double','double','double','double','double',...
                'string','string','string'},'VariableNames',{'width','height',...
                'slices','channels','frames','images','xRes','yRes','zRes',...
                'dataType','dimOrder','cOrder'});
            mopt_.dataType = "uint16";
            mopt_.dimOrder = string.empty(1,0);
            mopt_.cOrder = string.empty(1,0);

            mov = regmov(uint16.empty(), mopt_, [], []);
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
