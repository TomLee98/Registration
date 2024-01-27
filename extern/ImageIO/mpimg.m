classdef mpimg < handle
    %MEMMAPPER This is memmepfile interface defination
    % which is considered as value class, as same as matfile class but more
    % fast when mapped file is not big enough(<20GB)

    properties(Constant, Hidden)
        BYTES_UINT8 = 1
        BYTES_UINT16 = 2
        BYTES_UINT32 = 4
        BYTES_UINT64 = 8
    end

    properties(Access = private, Hidden)
        memptr          % 1-by-1 memmapfile object
        fmt             % 1-by-3 cell array, with {datatype, size, field}
        dimorder        % 1-by-d string array, dimension order of data
        isdisplay       % 1-by-1 logical, display flag
        isautoclear     % 1-by-1 logical, auto clear disk file flag
        isconst         % 1-by-1 logical, mark the file is constant or variable
        datatype        % 1-by-n char array, indicate data type, such as 'uint16', etc
    end

    properties(Access=public, Dependent)
        Data            % variable, get/set, (can't set if isconst is true)
        DimOrder        % variable, get/set, (can't set if isconst is true)
        IsAutoClear     % constant, get
        IsConst         % constant, get
        DataType        % variable, get
        DataSize        % variable, get
        DataDims        % variable, get
        DataBytes       % variable, get
    end

    methods
        function this = mpimg(folder_, file_, data_, dimorder_, autoclear_, const_, display_)
            %MEMMAPPER A Constructor
            % case 1: folder ~= [], file_ == [], data_ ~= []    % new file with random name
            % case 2: folder ~, file_ ~= [], data ~= []         % new file with given name
            % otherwise: throw exception: invalid use
            arguments
                folder_
                file_
                data_
                dimorder_  (1,:)    string = ["X","Y","C","Z","T"]
                autoclear_ (1,1)    logical = true
                const_     (1,1)    logical = false
                display_   (1,1)    logical = false     % for debugging
            end

            this.dimorder = dimorder_;
            this.isdisplay = display_;
            this.isautoclear = autoclear_;
            this.isconst = const_;

            if (~isempty(folder_) && isfolder(folder_)) ...
                    && isempty(file_) && ~isempty(data_)
                file_ = mpimg.genfilename(folder_);
                this.newmapfile(file_, data_);
            elseif isempty(folder_) && (~isempty(file_)&&~isfile(file_)) ...
                    && ~isempty(data_)
                this.newmapfile(file_, data_);
            else
                throw(MException("mpimg:invalidUse", ...
                    "Can not parse input arguments."));
            end
        end

        function r = get.Data(this)
            r = this.memptr.Data.mov;
        end

        function set.Data(this, r)
            if this.isconst == false
                if all(size(this.memptr.Data.mov) ...
                        == size(r))
                    % only change the value, write to disk
                    this.memptr.Data.mov = r;
                else
                    % remapping a new file, but keep the name
                    this.remaptr(r);
                end
            else
                throw(MException("mpimg:tryToModifyConstant", ...
                    "Constant value can not be modified."));
            end
        end

        function r = get.DimOrder(this)
            r = this.dimorder;
        end

        function set.DimOrder(this, r)
            arguments
                this
                r   (1,:) string = ["X","Y","C","Z","T"]
            end

            if numel(r) ~= numel(this.fmt)
                throw(MException("mpimg:invalidDimensions", ...
                    "Permute data need compelete dimension order."));
            end

            % permute the data
            [~, p] = ismember(r, this.fmt);

            if any(p==0)
                throw(MException("mpimg:invalidDimensionIndicator", ...
                    sprintf("Dimensions indicator must be in [%s]", this.fmt.join(","))));
            end
            if numel(unique(p)) ~= numel(p)
                throw(MException("mpimg:invalidDimensionIndicator", ...
                    sprintf("Permute data need compelete dimension order.")));
            end

            this.Data = permute(this.Data, p);  % easy calling

            this.dimorder = r;
        end

        function r = get.IsAutoClear(this)
            r = this.isautoclear;
        end

        function r = get.IsConst(this)
            r = this.isconst;
        end

        function r = get.DataType(this)
            r = this.datatype;
        end

        function r = get.DataSize(this)
            r = size(this.memptr.Data.mov);
        end

        function r = get.DataDims(this)
            r = ndims(this.memptr.Data.mov);
        end

        function r = get.DataBytes(this)
            switch this.datatype
                case {'uint8', 'int8'}
                    bytes_per_elem = this.BYTES_UINT8;
                case {'uint16','int16'}
                    bytes_per_elem = this.BYTES_UINT16;
                case {'uint32','int32','single'}
                    bytes_per_elem = this.BYTES_UINT32;
                case {'uint64','int64','double'}
                    bytes_per_elem = this.BYTES_UINT64;
            end
            r = prod(this.DataSize)*bytes_per_elem;
        end

        function delete(this)
            free(this);

            % calling default destructor
            % ~
        end

        function free(this)
            % this function will free the linked data
            % if autoclear is true, the disk file will be removed
            if ~isempty(this.memptr)
                file_ = this.memptr.Filename;

                % auto free by memmapfile->cleanup
                this.memptr = [];

                if this.isautoclear == true
                    % free disk file
                    delete(file_);
                end
            end
        end

        function link(this, file_, format_, dimorder_)
            % this function will link to an exist file
            arguments
                this
                file_
                format_    (1,:)    cell = cell(1,3)
                dimorder_  (1,:)    string = ["X","Y","C","Z","T"]
            end

            if (numel(format_) ~= 3) || ~ismember(format_{1}, ...
                    ["uint8","uint16","uint32","uint64","int8","int16", ...
                    "int32","int64","single","double"]) ...
                    || ~(isvector(format_{2})&&isPositiveIntegerValuedNumeric(format_{2})) ...
                    || ~isvarname(format_{3})
                throw(MException("mpimg:invalidFormatArg", ...
                    "Invalid file format argument."));
            end

            if isfile(file_)
                [~, ~, ext] = fileparts(file_);
                if ~strcmp(ext, '.dat')
                    throw(MException("mpimg:genmptr:invalidUse", ...
                        "Invalid file format."));
                end

                this.fmt = format_;
                this.dimorder = dimorder_;

                this.memptr = memmapfile(file_, ...
                    "Format", this.fmt, ...
                    "Writable", ~this.isconst);

                if this.isdisplay, fprintf("Memery mapping linked.\n"); end
            else
                throw(MException("mpimg:genmptr:invalidUse", ...
                    "File is not found or access permission denied."));
            end
        end

        function r = isempty(this)
            r = isempty(this.memptr);
        end
    end

    methods(Access = private)

        function remaptr(this, data_)
            % remove file and release resource
            file_ = this.memptr.Filename;

            if isfile(file_)
                if this.isautoclear
                    this.memptr = [];
                    delete(file_);
                else
                    throw(MException("mpimg:invalidFileName", ...
                        "Can not create a file has the same name with " + ...
                        "an already existed file when auto clear is " + ...
                        "not working."));
                end
            end

            % memmapping
            this.newmapfile(file_, data_);
        end

        function newmapfile(this, file_, data_)
            % memmapping
            msize = size(data_);
            this.datatype = class(data_);

            % replace the format
            this.fmt = {this.datatype, msize, 'mov'};

            if this.isdisplay, fprintf("Memery mapping -> %s ...\n", file_); end

            try
                fid = fopen(file_, "w");
                fwrite(fid, data_, this.datatype);
                fclose(fid);
            catch ME
                throwAsCaller(ME);
            end

            this.memptr = memmapfile(file_, ...
                "Format", this.fmt, ...
                "Writable", ~this.isconst);

            if this.isdisplay, fprintf("Memery mapping done.\n"); end
        end

    end

    methods(Static)
        function file_ = genfilename(folder_)
            % generate file name code: 18 chars
            code_idx = [randi(26,1,6)+64, randi(26,1,6)+96, randi(10,1,6)+47];
            code_idx = code_idx(randperm(18));

            % random suffix for avoiding same filename
            % and hidden the file information
            filename_ = char(code_idx);
            file_ = [folder_.char(), filesep, filename_, '.dat'];
        end
    end
end
