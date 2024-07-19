classdef mpimg < matlab.mixin.Copyable
    %MEMMAPPER This is memmepfile interface definition
    % which is considered as handle class, as same as mat-file class but more
    % fast when mapped file is not big enough(<20GB)
    % TODO: increment IO for low memory using 

    properties(Constant, Hidden)
        BYTES_UINT8 = 1
        BYTES_UINT16 = 2
        BYTES_UINT32 = 4
        BYTES_UINT64 = 8
        INNER_DIMENSION_ORDER = ["X","Y","C","Z","T"]
        INNER_BLOCK_SIZE = 50
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
        FileName        % variable, get
    end

    properties(Access = private, Hidden)
        fmt             % 1-by-3 cell array, with {datatype, size, field}
        dimorder        % 1-by-d string array, dimension order of data
        isdisplay       % 1-by-1 logical, display flag
        isautoclear     % 1-by-1 logical, auto clear disk file flag
        isconst         % 1-by-1 logical, mark the file is constant or variable
    end

    properties(Access = private, Hidden, NonCopyable)
        memptr          % 1-by-1 memmapfile object
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

            if any(~ismember(upper(dimorder_), this.INNER_DIMENSION_ORDER))
                throw(MException("mpimg:invalidDimensionIndicator", ...
                    "Unknown dimension indicator."));
            end
            if numel(unique(dimorder_)) ~= numel(dimorder_)
                throw(MException("mpimg:invalidDimensionIndicator", ...
                    "Dimension indicators must be unique."));
            end

            this.dimorder = upper(dimorder_);
            this.isdisplay = display_;
            this.isautoclear = autoclear_;
            this.isconst = const_;

            if (~isempty(folder_) && isfolder(folder_)) ...
                    && isempty(file_) && ~isempty(data_)
                file_ = mpimg.genfilename(folder_);
                this.newmptr(file_, data_);
            elseif isempty(folder_) && (~isempty(file_)&&~isfile(file_)) ...
                    && ~isempty(data_)
                this.newmptr(file_, data_);
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
                        == size(r)) ...
                        && strcmp(this.DataType, class(r))
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

            if numel(r) ~= numel(this.dimorder)
                throw(MException("mpimg:invalidDimensions", ...
                    "Permute data need compelete dimension order."));
            end

            % permute the data
            [~, p] = ismember(upper(r), this.dimorder);

            if any(p==0)
                throw(MException("mpimg:invalidDimensionIndicator", ...
                    sprintf("Dimensions indicator must be in [%s]", this.dimorder.join(","))));
            end
            if numel(unique(p)) ~= numel(p)
                throw(MException("mpimg:invalidDimensionIndicator", ...
                    sprintf("Permute data need compelete dimension order.")));
            end

            this.Data = permute(this.Data, p);  % easy calling

            this.dimorder = upper(r);
        end

        function r = get.IsAutoClear(this)
            r = this.isautoclear;
        end

        function r = get.IsConst(this)
            r = this.isconst;
        end

        function r = get.DataType(this)
            r = class(this.memptr.Data.mov);
        end

        function r = get.DataSize(this)
            r = size(this.memptr.Data.mov);
        end

        function r = get.DataDims(this)
            r = ndims(this.memptr.Data.mov);
        end

        function r = get.DataBytes(this)
            switch this.DataType
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

        function r = get.FileName(this)
            r = this.memptr.Filename;
        end
    end

    methods(Access = public)
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
                throw(MException("mpimg:link:invalidFormatArg", ...
                    "Invalid file format argument."));
            end

            if isfile(file_)
                this.linkextf(file_, format_, dimorder_);

                if this.isdisplay, fprintf("Memery mapping linked.\n"); end
            else
                throw(MException("mpimg:genmptr:invalidUse", ...
                    "File is not found or access permission denied."));
            end
        end

        function mptr = vcrop(this, dim_, r_)
            arguments
                this
                dim_  (1,1) string
                r_    (:,2) double {mustBePositive, mustBeInteger}
            end

            dim_ = upper(dim_).split("");
            dim_ = dim_(2:end-1);    % remove "" at front and tail
            [s, ~] = ismember(dim_, this.dimorder);
            if any(~s)
                throw(MException("mpimg:invalidCrop", ...
                    "No dimension for this operation."));
            end
            if numel(dim_) ~= size(r_, 1)
                throw(MException("mpimg:invalidCrop", ...
                    "Crop dimension not match."));
            end
            validateattributes(r_', "double", "nondecreasing");

            sz = repmat({':'}, 1, numel(this.INNER_DIMENSION_ORDER));
            [~, p] = ismember(dim_, this.INNER_DIMENSION_ORDER);
            for n = 1:numel(p)
                if r_(n, 2) > this.DataSize(dim_(n)==this.dimorder)
                    throw(MException("mpimg:invalidCropRange", ...
                        "Crop range is out of data size."));
                end
                sz{p(n)} = r_(n, :);
            end

            mptr = this.crop_(sz);
        end

        function mptr = tcrop(this, r_)
            arguments
                this
                r_    (1,:) double {mustBePositive, mustBeInteger}
            end

            % also call crop_
            if ~ismember("T", this.dimorder)
                throw(MException("mpimg:invalidCrop", ...
                    "No dimension for this operation."));
            end

            validateattributes(r_, "double", "nondecreasing");

            sz = repmat({':'}, 1, numel(this.INNER_DIMENSION_ORDER));
            [~, p] = ismember("T", this.INNER_DIMENSION_ORDER);
            if r_(end) > this.DataSize("T"==this.dimorder)
                throw(MException("mpimg:invalidCropRange", ...
                    "Crop range is out of data size."));
            end

            % transform to char array representation
            r_ = "[" + string(r_).join(",") + "]";
            r_ = char(r_);

            sz{p} = r_;

            mptr = this.crop_(sz);
        end

        function promote(this, dim_)
            % This function promotes a new dimension next the last dimension
            % for example, [X,Y,Z] could be promoted as [X,Y,Z,T], etc
            arguments
                this
                dim_    (1,1)   string
            end

            if ~ismember(dim_, this.INNER_DIMENSION_ORDER)
                throw(MException("mpimg:promote:invalidPromote", ...
                    "Only 'X','Y','C','Z','T' are valid "))
            end

            % TODO: 
        end

        function ctset(this, v, cr, tr)
            % This function selects cr and tr at C & T dimension, then
            % modify the data by v
            arguments
                this
                v   (:,:,:,:)        % dim as [x,y,z,t]
                cr  (1,:)   logical
                tr  (1,:)   double  {mustBePositive, mustBeInteger}
            end

            % optimize the usual case for fast saving
            switch this.dimorder.join("")
                case "XYCZT"
                    this.memptr.Data.mov(:,:,cr,:,tr) = v;
                case "XYZCT"
                    this.memptr.Data.mov(:,:,:,cr,tr) = v;
                otherwise
                    % generate the crop expression
                    expr = this.gen_expstr(cr, tr, "v");

                    % simple modifying
                    eval(expr);
            end
        end

        function fliplr(this)
            % This function flips left and right on XY plane at the same
            % place
            this.flip_("lr");
        end

        function flipud(this)
            % This function flips up and down on XY plane at the same 
            % place
            this.flip_("ud");
        end

        function flipbt(this)
            % This function flips up and down on XZ plane at the same 
            % place
            this.flip_("bt");
        end

        function rot90(this, k)
            % This function rotate XY plane 90 degrees k times
            this.rot90_(k);
        end

        function this = plus(this, r_)
            % operator '+' overloading, data plus number on each channel

            arguments
                this
                r_  double  {mustBeVector}
            end

            nc = this.DataSize("C"==this.DimOrder);

            if numel(r_) ~= nc
                throw(MException("regmov:plus:invalidChannelsNumber", ...
                    "Channels number not match."));
            end
            
            % range splitter for blocked data saving
            nt = this.DataSize("T"==this.DimOrder);
            t_range_ = range_splitter(sprintf("1:%d", nt), this.INNER_BLOCK_SIZE);

            sz_ = repmat({':'}, 1, numel(this.INNER_DIMENSION_ORDER));

            for c = 1:nc
                sz_{"C"==this.dimorder} = c;

                for t = 1:numel(t_range_)
                    sz_{"T"==this.dimorder} = t_range_(t);

                    expr = "this.memptr.Data.mov("+string(sz_).join(",")+")";
                    expr = sprintf("%s=%s+%.1f;", expr, expr, r_(c));

                    eval(expr);
                end
            end
        end

        function this = minus(this, r_)
            % operator '-' overloading, data plus number on each channel
            arguments
                this
                r_  double  {mustBeVector}
            end

            this = this.plus(-r_);
        end

        function r = isempty(this)
            r = isempty(this.memptr);
        end

        % destructor
        function delete(this)
            free(this);

            % calling default destructor
            % ~
        end
    end

    methods(Access = protected)
        % Override copyElement method:
        function cpt = copyElement(this)
            % Make a shallow copy of all properties except memptr
            cpt = copyElement@matlab.mixin.Copyable(this);

            % Make a deep copy of the memptr object
            [folder_, ~, ~] = fileparts(this.memptr.Filename);
            file_ = mpimg.genfilename(folder_);
            try
                fid = fopen(file_, "w");
                fwrite(fid, this.Data, this.DataType);
                fclose(fid);
            catch ME
                throwAsCaller(ME);
            end
            cpt.memptr = memmapfile(file_, ...
                "Format", this.fmt, ...
                "Writable", ~this.isconst);
        end
    end

    methods(Access = private, Hidden)

        % This function remaps this to a new file which contains data_
        % and keeps this io properties and file name
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
            this.newmptr(file_, data_);
        end

        % This function create a new file contains data_
        function newmptr(this, file_, data_)
            % memmapping
            msize_ = size(data_);
            datatype_ = class(data_);

            % replace the format
            this.fmt = {datatype_, msize_, 'mov'};

            if this.isdisplay, fprintf("Memery mapping -> %s ...\n", file_); end

            try
                fid = fopen(file_, "w");
                fwrite(fid, data_, datatype_);
                fclose(fid);
            catch ME
                throwAsCaller(ME);
            end

            this.memptr = memmapfile(file_, ...
                "Format", this.fmt, ...
                "Writable", ~this.isconst);

            if this.isdisplay, fprintf("Memery mapping done.\n"); end
        end

        function linkextf(this, file_, format_, dimorder_)
            [~, ~, ext] = fileparts(file_);
            if ~strcmp(ext, '.dat')
                throw(MException("mpimg:linkextmptr:invalidUse", ...
                    "Can not link to invalid file with suffix is not .dat."));
            end

            this.fmt = format_;
            this.dimorder = dimorder_;

            this.memptr = memmapfile(file_, ...
                "Format", this.fmt, ...
                "Writable", ~this.isconst);
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

        % This function crop data and return another mpimg object, which
        % contains the cropped data
        function mptr = crop_(this, sz_)
            % sz must has the same dimension with dimorder
            if numel(sz_) ~= numel(this.INNER_DIMENSION_ORDER)
                throw(MException("mpimg:crop:invalidUse", ...
                    "Dimension not match."));
            end

            % parse the cropped range
            for n = 1:numel(sz_)
                if ~ischar(sz_{n})
                    sz_{n} = [num2str(sz_{n}(1)),':',num2str(sz_{n}(2))];
                end
            end
            [~, p] = ismember(this.dimorder, this.INNER_DIMENSION_ORDER);
            sz_ = sz_(p);

            % fill time span automatically if ':' for compatibility
            if strcmp(sz_{end}, ':')
                sz_{end} = ['1',':', string(this.DataSize("T"==this.DimOrder)).char()];
            end

            % range splitter for blocked data saving
            t_range_ = range_splitter(sz_{end}, this.INNER_BLOCK_SIZE);

            % create the primer
            sz_{end} = t_range_(1);
            expr = "this.Data("+string(sz_).join(",")+");";
            data_ = eval(expr);
            % generate new memmapping
            [folder_, ~, ~] = fileparts(this.memptr.Filename);
            mptr = mpimg(string(folder_), [], data_, ...
                this.dimorder, this.isautoclear, this.isconst, this.isdisplay);

            % link the data blocks to primer
            for k = 2:numel(t_range_)
                sz_{end} = t_range_(k);
                expr = "this.Data("+string(sz_).join(",")+");";
                data_ = eval(expr);

                mptr.concat(data_);
            end

            if this.isdisplay == true
                fprintf("New mpimg object inherits this properites.")
            end
        end

        function flip_(this, direct)
            dims_T_loc = find(this.dimorder=="T");
            NT = this.DataSize(dims_T_loc);
            Group_N = ceil(NT/this.INNER_BLOCK_SIZE);

            [~, pdim_xy] = ismember(this.dimorder, this.INNER_DIMENSION_ORDER);
            [~, pdim_xz] = ismember(this.dimorder, ["X","Z","C","Y","T"]);
            [~, invpdim_xy] = ismember(this.INNER_DIMENSION_ORDER, this.dimorder);
            [~, invpdim_xz] = ismember(["X","Z","C","Y","T"], this.dimorder);

            % use subsasgn structure
            expr = repmat({':'}, 1, this.DataDims);

            for k = 1:Group_N
                % read the data along dimension T
                t_r = [(k-1)*this.INNER_BLOCK_SIZE+1, min(k*this.INNER_BLOCK_SIZE, NT)];
                expr{dims_T_loc} = string(t_r).join(":");
                expstr = sprintf("this.memptr.Data.mov(%s)", string(expr).join(","));
                D = eval(expstr);

                % permute D and flip
                switch direct
                    case "lr"
                        D = fliplr(permute(D, pdim_xy));   % fliplr at dim: XYCZT

                        % permute reverse
                        D = permute(D, invpdim_xy);        %#ok<NASGU>
                    case "ud"
                        D = flipud(permute(D, pdim_xy));   % flipud at dim: XYCZT

                        % permute reverse
                        D = permute(D, invpdim_xy);        %#ok<NASGU>
                    case "bt"
                        D = flipud(permute(D, pdim_xz));   % flipud at dim: XZCYT

                        % permute reverse
                        D = permute(D, invpdim_xz);        %#ok<NASGU>
                    otherwise
                end
                
                % save the data along dimension T
                eval(expstr+" = D;");
            end
        end

        function rot90_(this, N)
            % find temporal dimension and generate extraction prompt
            sz_ = repmat({':'}, 1, this.DataDims);
            locT = ("T"==this.DimOrder);
            % TODO: locT matches this.Data dimension order?
            sz_{locT} = ['1',':', string(this.DataSize(locT)).char()];

            % range splitter for blocked data saving
            t_range_ = range_splitter(sz_{locT}, this.INNER_BLOCK_SIZE);

            % create the primer
            sz_{locT} = t_range_(1);
            expr = "this.Data("+string(sz_).join(",")+");";
            data_ = eval(expr);

            % rotation
            data_ = rot90(data_, N);

            % generate new memmapping
            [folder_, ~, ~] = fileparts(this.memptr.Filename);
            mptr = mpimg(string(folder_), [], data_, ...
                this.dimorder, this.isautoclear, this.isconst, this.isdisplay);

            % link the data blocks to primer
            for k = 2:numel(t_range_)
                sz_{locT} = t_range_(k);
                expr = "this.Data("+string(sz_).join(",")+");";
                data_ = eval(expr);
                data_ = rot90(data_, N);

                mptr.concat(data_);
            end

            if this.isdisplay == true
                fprintf("New mpimg object inherits this properites.")
            end
        end

        function concat(this, data_)
            % this function concatenate the data on the last dimension,
            % which requires the compatible data size
            % Note that: concatenate is self operation for avoiding IO overhead

            sz_ = size(data_);

            % validate dimensionality
            if numel(sz_) > numel(this.dimorder) ...
                    || numel(sz_) < numel(this.dimorder) - 1
                throw(MException("mpimg:concatnate:invalidUse", ...
                    "Data dimensionality is not compatible."));
            end

            % validate size compatibility
            if (numel(sz_)==numel(this.dimorder) && any(sz_(1:end-1)~=this.DataSize(1:end-1))) ...
                    || (numel(sz_)==numel(this.dimorder)-1 && any(sz_~=this.DataSize(1:end-1)))
                throw(MException("mpimg:concatnate:invalidUse", ...
                    "Data size is not compatible."));
            end

            if class(data_) ~= this.DataType
                warning("mpimg:concatnate:badUse", ...
                    "Data format will be changed to %s.", this.DataType);
                data_ = cast(data_, this.DataType);
            end

            % record this file name
            file_ = this.FileName;
            datatype_ = this.DataType;
            fmt_ = this.DataSize;

            % cut the memptr linkage
            this.memptr = [];

            % open the file with 'append' mode and write
            try
                fid = fopen(file_, "a");
                fwrite(fid, data_, datatype_);
                fclose(fid);
            catch ME
                throwAsCaller(ME);
            end

            % update and re-link the new file
            if numel(sz_) == numel(this.dimorder) - 1
                nb = 1;
            else
                nb = sz_(end);
            end
            msize_ = [fmt_(1:end-1), fmt_(end)+nb];

            format_ = {datatype_, msize_, 'mov'};

            this.linkextf(file_, format_, this.dimorder);
        end

        function expr = gen_expstr(this, cr, tr, istr)
            t_loc = find(this.dimorder=="T");
            c_loc = find(this.dimorder=="C");
            cr = string(find(cr));

            if isempty(cr) || isempty(t_loc) || isempty(c_loc)
                throw(MException("regmov:ctset:invalidMovie", ...
                    "Bad calling: invalid movie dimension."));
            end

            expr = "";

            for dp = 1:this.DataDims
                if dp == c_loc
                    expr = expr + cr;
                elseif dp == t_loc
                    expr = expr + tr;
                else
                    expr = expr + ":";
                end
                if dp ~= this.DataDims, expr = expr + ","; end
            end

            expr = sprintf("this.memptr.Data.mov(%s)=%s;", expr, istr);
        end
    end

    methods(Static, Hidden)
        function file_ = genfilename(folder_)
            % generate file name code: 18 chars
            code_idx = [randi(26,1,6)+64, randi(26,1,6)+96, randi(10,1,6)+47];
            code_idx = code_idx(randperm(18));

            % random suffix for avoiding same filename
            % and hidden the file information
            filename_ = char(code_idx);

            if isstring(folder_), folder_ = folder_.char(); end

            file_ = [folder_, filesep, filename_, '.dat'];
        end
    end

    methods(Static)
        function folder_ = findtmpfolder(volopt_, numhistory_)
            %FINDTMPFOLDER This function find a temporary folder in a disk part, which
            %could contains the memory mapping files
            % input:
            %   - volopt: 1-by-12 table, with image options
            %   - numhistory: 1-by-1 positive integer, the number of operation history
            %                   list, for space estimation
            % output:
            %   - tf: 1-by-1 string, the temporary folder
            arguments
                volopt_      (1,12)  table
                numhistory_  (1,1)   double {mustBeNonnegative, mustBeInteger} = 0
            end

            import java.io.*;
            import javax.swing.filechooser.*;
            import java.lang.*;

            warning('off', 'MATLAB:MKDIR:DirectoryExists');

            if ispc()
                files_ = File.listRoots();
                user_name = string(System.getProperty("user.name"));
                disk_table = table('Size',[numel(files_), 5], ...
                    'VariableTypes',{'string','double','logical','logical','logical'}, ...
                    'VariableNames',{'letter','free_size','readable','writable','is_local'});
                for p = 1:numel(files_)
                    disk_table.letter(p) = string(files_(p).getPath());
                    disk_table.free_size(p) = files_(p).getFreeSpace();
                    disk_table.readable(p) = files_(p).canRead();
                    disk_table.writable(p) = files_(p).canWrite();
                    disk_table.is_local(p) = is_local_disk(files_(p));
                end

                % find the maximum disk part as temporary root
                vp = disk_table.readable & disk_table.writable ...
                    & disk_table.is_local;
                disk_table = disk_table(vp, :);
                [fs, idx] = max(disk_table.free_size);
                if fs <= calc_volsize(volopt_, numhistory_)
                    throw(MException("findtmpfolder:noValidSpace", ...
                        "Data is too big, try to reduce number of histories or using " + ...
                        "mpimgs object instead."));
                end

                folder_ = [disk_table.letter(idx).char(), 'Reg3DCache'];
               
                try
                    mkdir(folder_);
                    fileattrib(folder_, '+h', '', 's');
                    folder_ = [folder_, filesep, user_name.char()];
                    mkdir(folder_);
                catch ME
                    switch ME.identifier
                        case "matlab.io.FileError"
                            % ??
                    end
                end

            elseif isunix()
                if isfolder('/data/.Reg3DCache')
                    % create user folder on /data/Reg3DCache
                    user_name = string(System.getProperty("user.name"));
                    folder_ = ['/data/.Reg3DCache/', char(user_name)];
                    try
                        mkdir(folder_);
                    catch ME
                        throwAsCaller(ME);
                    end
                else
                    throw(MException("Cache folder lost. " + ...
                        "Please connect to the administrator."));
                end
            end

            warning('on', 'MATLAB:MKDIR:DirectoryExists');

            function sz_ = calc_volsize(volopt_, numhistory_)
                switch volopt_.dataType
                    case {'uint8', 'int8'}
                        bytes_per_elem = 1;
                    case {'uint16','int16'}
                        bytes_per_elem = 2;
                    case {'uint32','int32','single'}
                        bytes_per_elem = 4;
                    case {'uint64','int64','double'}
                        bytes_per_elem = 8;
                end
                sz_ = volopt_.width * volopt_.height * volopt_.images ...
                    * bytes_per_elem * (numhistory_+1);
            end

            function tf_ = is_local_disk(file_)
                % NOTE: This function can't distinguish movable disk and 
                % local disk

                % construct java object
                lang = System.getProperty("user.language");
                fileSystemView = FileSystemView.getFileSystemView();

                if string(lang) == "zh"
                    marker = "本地磁盘";
                else 
                    marker = "LocalDisk"; % ?? not confirm
                end
                diskType = fileSystemView.getSystemTypeDescription(file_);
                tf_ = (string(diskType) == marker);
            end
        end
    end
end

% ===================== utility function ==========================
function nbs = range_splitter(s, bsz)
% This function help to split the range string as number-block setting
% Input:
%   - s: the sliced indices indicator string, like "1:20" or "1,3,5,7,8" etc
%   - bsz: 1-by-1 positive integer, the splitted block size
% Output:
%   - nbs: q-by-1 string array, with splitted indicator s

v = str2num(s); %#ok<ST2NM>

nb = ceil(numel(v) / bsz);
nbs = strings(nb, 1);
nm = numel(v);

% soop: string object oriented programming(
for k = 1:nb
    nbs(k) = string(v((k-1)*bsz+1:min(k*bsz, nm))).join(",");
    nbs(k) = "[" + nbs(k) + "]";
end

end
