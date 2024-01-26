classdef mpimg < handle
    %MEMMAPPER This is memmepfile interface defination
    properties(Constant, Hidden)
        INNER_DIMENSION = ["X","Y","C","Z","T"]
        BYTES_UINT16 = 2
    end

    properties(Access = private, Hidden)
        mptr
        display
    end

    properties(GetAccess=public, Dependent)
        Filename
        Data
        Size
        Bytes
    end

    properties(GetAccess=public, SetAccess=private)
        Dim
    end
    
    methods
        function this = mpimg(folder_, data_, dim_, display_)
            %MEMMAPPER A Constructor
            arguments
                folder_  (1,1)        string {mustBeFolder}
                data_    (:,:,:,:,:)  uint16
                dim_     (1,5)        string = ["X","Y","C","Z","T"]
                display_ (1,1)        logical = true
            end

            this.Dim = dim_;
            this.display = display_;

            this.genmptr(folder_, [], data_);
        end

        function r = get.Filename(this)
            r = this.mptr.Filename;
        end

        function r = get.Data(this)
            fn = this.mptr.Format{3};
            r = this.mptr.Data.(fn);
            [~, p] = ismember(this.Dim, this.INNER_DIMENSION);
            r = permute(r, p);
        end

        function r = get.Size(this)
            [~, p] = ismember(this.Dim, this.INNER_DIMENSION);
            r = this.mptr.Format{2}(p);
        end

        function r = get.Bytes(this)
            r = prod(this.Size)*this.BYTES_UINT16;
        end
        
        function this = SetData(this, data_, dim_)
            %SetData This function delete the old file and generate a new
            %file but with the same name
            arguments
                this
                data_  (:,:,:,:,:)  uint16
                dim_   (1,5)        string = ["X","Y","C","Z","T"]
            end

            % 1. get the filename
            filename = this.Filename;

            % 2. remove mptr for resource release
            this.mptr = [];

            % 3. delete the old file
            delete(filename);

            % 4. create new pointer
            this.Dim = dim_;

            display_ = this.display;
            this.display = false;

            this.genmptr([], filename, data_);

            this.display = display_;
        end

        function r = GetSliceAt(this, cidx, zidx, vidx)
            r = squeeze(this.Data(:,:,cidx, zidx,vidx));
        end

        function r = GetVolumeAt(this, cidx, vidx)
            r = squeeze(this.Data(:,:,cidx,:,vidx));
        end

        function r = GetChannelAt(this, cidx)
            r = squeeze(this.Data(:,:,cidx,:,:));
        end

        function r = GetMptr(this)
            r = this.mptr;
        end

        % =================== operator overloading ==================

        function delete(this)
            filename = this.Filename;

            % remove mptr for resource release
            this.mptr = [];

            % free disk file
            delete(filename);
        end

        function ind = end(this, k, n)
            if k < n
                ind = this.Size(k);
            else
                ind = prod(this.Size(k:end));
            end
        end

        function r = subsref(this, s)
            arguments
                this
                s   (1,1)   struct
            end

            % inner class callback
            % TODO: dot calls member function with data
            switch s.type
                case '()'
                    % do crop on data
                    r = subsref(this.Data, s);
                case '.'
                    r = this.(s.subs);
                otherwise
                    throw(MException("Memmapper:invalidOperation", ...
                        "Not overloaded operator."));
            end
        end

        function this = subsasgn(this, s, b)
            arguments
                this
                s   (1,1) struct
                b
            end

            switch s.type
                case '()'
                    r = subsasgn(this.Data, s, b);
                    fn = this.mptr.Format{3};
                    this.mptr.Data.(fn) = r;
                otherwise
                    throw(MException("Memmapper:invalidOperation", ...
                        "Dot operator can't set left value."));
            end
        end

        function r = ndims(this)
            r = numel(this.Size);
        end

        function r = size(this, dims)
            arguments
                this
                dims double {mustBePositive, mustBeInteger} ...
                    = 1:numel(this.Size)
            end

            r = this.Size(dims);
        end

        function this = plus(this, d)
            fn = this.mptr.Format{3};
            % write to disk
            this.mptr.Data.(fn) ...
                = this.mptr.Data.(fn) + d;
        end

        function this = minus(this, d)
            fn = this.mptr.Format{3};
            % write to disk
            this.mptr.Data.(fn) ...
                = this.mptr.Data.(fn) - d;
        end
    end

    methods(Access = private)
        function this = genmptr(this, folder_, file_, data_)
            % permute the dimension
            if any(this.Dim ~= this.INNER_DIMENSION)
                [~, p] = ismember(this.INNER_DIMENSION, this.Dim);
                data_ = permute(data_, p);
            end

            msize = size(data_);

            if isempty(file_)
                % save mov and get the file pointer
                % generate file name code: 18 chars
                code_idx = [randi(26,1,6)+64, randi(26,1,6)+96, randi(10,1,6)+47];
                code_idx = code_idx(randperm(18));

                % random suffix for avoiding same filename
                % and hidden the file information
                filesuffix = char(code_idx);
                file_ = [folder_.char(), filesep, filesuffix, '.dat'];
            end

            if this.display
                fprintf("Memery mapping -> %s ...\n", file_);
            end

            try
                fid = fopen(file_, "w");
                fwrite(fid, data_, "uint16");
                fclose(fid);
            catch ME
                switch ME.identifier
                    case "MATLAB:FileIO:InvalidFid"
                        throw(MException("Memmapper:badPermission", ...
                            "You may haven't the delete permission."));
                    otherwise
                        throwAsCaller(ME);
                end
            end

            if this.display
                fprintf("Memery mapping done.\n");
            end

            this.mptr = memmapfile(file_, ...
                "Format",{'uint16', msize, 'mov'}, ...
                "Writable", true);      % value modification
        end
    end
end
