classdef regmov
    %REGMOV This class is regmov defination, which construct the key data
    %struct for Register

    properties(Constant, Hidden)
        
    end

    properties(Access=private, Hidden)
        mptr            % 1-by-1 mpimg/mpimgs object
        mopt            % 1-by-12 table, with movie options
        t               % t-by-1 double, camera time
        zproj           % double, z-projection of mptr.Data
        zprojalg        % 1-by-1 string, method of z-projection
        tform           % t-by-3 cell, with {global, local, manual} transformation unit
        bkg             % c-by-1 double, c for color channel number
    end

    properties(GetAccess=public, Dependent)
        Movie           % variable, get/set, manual update
        MetaData        % variable, get/set, manual update
        Time            % variable, get/set, manual update
        MovieZProj      % variable, get    , auto update
        ZProjMethod     % variable, set/get, manual update
        Transformation  % variable, get/set, manual update
        Bytes           % variable, get    , /
        Background      % variable, get/set, manual update
    end

    
    methods
        function this = regmov(mptr_, mopt_, t_, bkg_)
            %REGMOV A constructor
            arguments
                mptr_   (1,1)
                mopt_   (1,12)  table  {mustBeMovieOptions}
                t_      (:,1)   double {mustBeNonnegative}
                bkg_    (1,:)   double {mustBeNonnegative}
            end

            if ~ismember(class(mptr_), ["mpimg", "mpimgs"])
                throw(MException("regmov:invalidMoviePointer", ...
                    "Invalid movie pointer, Only mpimg or mpimgs object is supported."));
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
            this.bkg = bkg_;
        end

        function r = get.Movie(this)
            % this operation will call deep copy (constructor)
            r = this.mptr;
        end

        function this = set.Movie(this, r_)
            arguments
                this
                r_  (1,1)
            end

            if ~ismember(class(r_), ["mpimg", "mpimgs"])
                throw(MException("regmov:invalidMoviePointer", ...
                    "Invalid movie pointer, Only mpimg or mpimgs object is supported."));
            end

            this.mptr = r_;
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
            [vz, pz] = ismember("Z", this.mptr.DimOrder);
            if ~vz
                throw(MException("regmov:invalidOperation", ...
                    "Movie does not have dimension: 'Z'."));
            end
            sz_noz = this.mptr.DataSize;
            sz_noz(pz) = [];
            
            if any(size(this.zproj)~=sz_noz)
                % update z projection automatically
                this.zproj = Projection(this.mptr.Data, this.zprojalg, ...
                    pz);
            end

            r = this.zproj;
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
            zproj_ = this.zproj; %#ok<NASGU>
            prop_ = whos("zproj_");
            r.mem = prop_.Bytes;
            r.map = this.mptr.DataBytes;
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

    methods(Static)
        function mustBeMovieOptions(A)
            arguments
                A   (1, 12) table
            end
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

