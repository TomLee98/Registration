classdef regmov
    %REGMOV This class is regmov defination, which construct the key data
    %struct for Register

    properties(Constant, Hidden)
        BYTES_UINT8 = 1
        BYTES_UINT16 = 2
        BYTES_UINT32 = 4
        BYTES_UINT64 = 8
    end

    properties(Access=private, Hidden)
        mptr            % 1-by-1 mpimg/mpimgs object or numeric array
        mopt            % 1-by-12 table, with movie options
        t               % t-by-1 double, camera time
        zprojalg        % 1-by-1 string, method of z-projection
        tform           % t-by-3 cell, with {global, local, manual} transformation unit
        bkg             % c-by-1 double, c for color channel number
    end

    properties(Access=public, Dependent)
        Movie           % variable, get/set, manual update
        MetaData        % variable, get/set, manual update
        Time            % variable, get/set, manual update
        MovieZProj      % variable, get    , /
        ZProjMethod     % variable, set/get, manual update
        Transformation  % variable, get/set, manual update
        Bytes           % variable, get    , /
        Background      % variable, get/set, manual update
    end

    
    methods
        function this = regmov(mptr_, mopt_, t_, bkg_)
            %REGMOV A constructor
            arguments
                mptr_   
                mopt_   (1,12)  table  {regmov.mustBeMovieOptions}
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
            this.tform = cell(0,3);
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

        function this = set.Movie(this, r_)
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

        function this = set.MetaData(this, r_)
            arguments
                this
                r_  (1,12)  table
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

        function this = set.Time(this, r_)
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

        function this = set.ZProjMethod(this, r_)
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

        function this = set.Transformation(this, r_)
            arguments
                this
                r_  (:, 3) cell
            end
            if numel(r_) ~= this.mopt.frames
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

            if ismember(class(mptr_), ["mpimg", "mpimgs"])
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

        function this = set.Background(this, r_)
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

    methods(Static, Hidden)
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

            validateattributes(A.width, "double", {'scalar','integer','positive'});
            validateattributes(A.height, "double", {'scalar','integer','positive'});
            validateattributes(A.slices, "double", {'scalar','integer','positive'});
            validateattributes(A.channels, "double", {'scalar','integer','positive'});
            validateattributes(A.frames, "double", {'scalar','integer','positive'});
            validateattributes(A.images, "double", {'scalar','integer','positive'});
            validateattributes(A.xRes, "double", {'scalar','positive'});
            validateattributes(A.yRes, "double", {'scalar','positive'});
            validateattributes(A.zRes, "double", {'scalar','positive'});
            validateattributes(A.dataType, "string", {'scalar'});
            validateattributes(A.dimOrder, "string", {'vector'});
            validateattributes(A.cOrder, "string", {'vector'});
        end
    end
end

