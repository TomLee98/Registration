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
        Background      % variable, get/set, stored
    end

    methods
        function this = regmov(mptr_, mopt_, t_, bkg_)
            %REGMOV A constructor
            arguments
                mptr_   
                mopt_   (1,12)  table  {mustBeMovieOptions}
                t_      (:,1)   double {mustBeNonnegative}
                bkg_    (1,:)   double {mustBeNonnegative}
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
                r_  (1,12)  table   {regmov.mustBeMovieOptions}
            end

            this.mopt = r_;
        end

        function r = get.Time(this)
            if numel(this.t) ~= this.mopt.frames
                throw(MException("regmov:badTime", ...
                    "Number of Timesteps not match."));
            end

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
                throw(MException("regmov:badTime", ...
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

        function set.Background(this, r_)
            arguments
                this
                r_  (1,:)   double {mustBeNonnegative}
            end
            if numel(r_) ~= numel(this.mopt.cOrder)
                throw(MException("regmov:badTime", ...
                    "Number of color channels not match."));
            end

            this.Background = r_;
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

    methods(Access=public, Hidden)
        function r = isempty(this)
            r = isempty(this.Movie);
        end

        function movobj = crop(this, dim_, r_)
            arguments
                this
                dim_  (1,1) string
                r_    (:,2) double {mustBePositive, mustBeInteger}
            end

            mopt_ = this.mopt;
            t_ = this.t;
            tf_ = this.tform;
            bkg_ = this.bkg;

            if ismember(class(this.mptr), ["mpimg", "mpimgs"])
                mptr_ = this.mptr.crop(dim_, r_);
                dim_ = dim_.split("");
                dim_ = dim_(2:end);
                for n = 1:numel(dim_)
                    switch dim_(n)
                        case "X"
                            mopt_.width = mopt_.width - (diff(r_(n,:))+1);
                        case "Y"
                            mopt_.height = mopt_.height - (diff(r_(n,:))+1);
                        case "C"
                            mopt_.channels = mopt_.channels - (diff(r_(n,:))+1);
                            bkg_ = bkg_(r_(n,1):r_(n:2));
                        case "Z"
                            mopt_.slices = mopt_.slices - (diff(r_(n,:))+1);
                        case "T"
                            mopt_.frames = mopt_.frames - (diff(r_(n,:))+1);
                            t_ = t_(r_(n,1):r_(n:2));
                            tf_ = tf_(r_(n,1):r_(n:2), :);
                    end
                end
                mopt_.images = mopt_.slices*mopt_.channels*mopt_.frames;

                movobj = regmov(mptr_, mopt_, t_, bkg_);
                movobj.Transformation = tf_;
            elseif isnumeric(this.mptr)
                if numel(dim_) ~= size(r_, 1)
                    throw(MException("regmov:invalidCrop", ...
                        "Crop dimension not match."));
                end
                validateattributes(r_, "double", "nondecreasing");
                %TODO: ANY DIMENSION ORDER SUPPORT
                warning("regmov:unfinishFunction", "Only partial functions " + ...
                    "are supported.");
                switch dim_
                    case "XY"
                        mov = CropXY(this.Movie, r_');
                        mopt_.width = mopt_.width - (diff(r_(1,:))+1);
                        mopt_.height = mopt_.height - (diff(r_(2,:))+1);
                    case "Z"
                        mov = CropZ(this.Movie, r_);
                        mopt_.slices = mopt_.slices - (diff(r_)+1);
                    case "T"
                        mov = CropT(this.Movie, r_);
                        mopt_.frames = mopt_.slices - (diff(r_)+1);
                        t_ = t_(r_(1):r_(2));
                        tf_ = tf_(r_(1):r_(2));
                    otherwise
                        throw(MException("regmov:invalidOperation", ...
                            "Unsupported operation."));
                end
                mopt_.images = mopt_.slices*mopt_.channels*mopt_.frames;

                movobj = regmov(mov, mopt_, t_, this.bkg);
                movobj.Transformation = tf_;
            end
        end

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
                this.mptr = mpimg(tmpfolder, [], this.mptr, this.mopt.dimOrder);
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