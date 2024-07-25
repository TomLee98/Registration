classdef ImageReader < handle
    %IMAGELOADER This class for easy data loading, which hides the low
    %level readers working, as an adapter between files reader and regmov
    %class

    properties(Constant, Hidden)
        INNER_DIM_ORDER = ["X","Y","C","Z","T"]
    end
    
    properties(GetAccess=public, Dependent)
        Data
        T
        MetaData
        FileName
        State
    end

    properties(SetAccess=immutable, Hidden)
        caller          % calling object, must has SetProgressBar method
        mapflag         % 1-by-1 logical, memmap pointer flag

    end

    properties(Access=private, Hidden)
        data
        t
        metadata
        folder      (1,1)   string = ""
        srcfile     (1,1)   string = ""
        state       (1,1)   string  {mustBeMember(state, ["on", "off"])} = "off"
    end

    methods
        function this = ImageReader(caller_, mapflag_)
            %IMAGELOADER A Constructor
            arguments
                caller_     (1,1)   Register
                mapflag_    (1,1)   logical = false
            end
            this.caller = caller_;
            this.mapflag = mapflag_;
        end

        function r = get.Data(this)
            r = this.data;
        end

        function r = get.T(this)
            r = this.t;
        end

        function r = get.MetaData(this)
            r = this.metadata;
        end

        function r = get.FileName(this)
            r = this.srcfile;
        end

        function r = get.State(this)
            r = this.state;
        end

    end

    methods(Access=public)
        function load(this, fname, block)
            % This function load the data by calling the image reader
            arguments
                this
                fname   (1,1)   string
                block   (1,1)   double {mustBePositive, mustBeInteger} = 50;
            end

            if ~isfile(fname)
                % uigetfile seletion
                [dstfile, path] = uigetfile(...
                    {'*.ims;*.nd2;*.tif','microscope volume files(*.ims,*.tif,*.nd2)'},...
                    'volumes series selector');
                if isequal(dstfile,0) || isequal(path,0)
                    this.state = "off";
                    return;
                else
                    this.srcfile = string(fullfile(path, dstfile));
                end
            else
                this.srcfile = fname;
            end

            [~, ~, ext] = fileparts(this.srcfile);
            % import the data loader setting
            pathToReader = fileparts(mfilename('fullpath'));
            rf = importdata([pathToReader, filesep,'configuration.ini']);
            rf = string(rf).split(":");
            if ~ismember(upper(ext), rf(:,1))
                throw(MException("ImageLoader:invalidImageLoadingFunc", ...
                        "Can not parse function without registration."));
            end

            % pfunc: [meta, data] = pfunc(file, tspan)
            pfunc = str2func(rf(upper(ext)==rf(:,1), 2));

            % if turbo loading enabled
            ft = ImageReader.fast_loading(ext);

            % load the metadata
            info = loadmeta_bf(this.srcfile);
            this.metadata = info.opts;
            this.t = info.rt;

            this.folder = mpimg.findtmpfolder(this.metadata, 0);

            switch lower(ext)
                case ".tif"
                    if ft == true
                        % reset the block size
                        block = this.metadata.frames;
                        % indeterminate progress bar
                        this.caller.SetProgressBar(0, true);
                    end
                otherwise
                    block = min(ceil(this.metadata.frames/20), block);
                    this.caller.SetProgressBar(0);
            end

            % get data piecewiselly
            n_piece = ceil(this.metadata.frames/block);

            % turn on loader
            this.state = "on";

            if this.mapflag == true
                % generate a new temporary file
                dstfile = mpimg.genfilename(this.folder);
                fid = fopen(dstfile, "a");
                for k = 1:n_piece
                    if this.state == "off"
                        break;
                    end
                    tspan = [(k-1)*block+1, min(k*block, this.metadata.frames)];
                    % load piecewise data
                    img = pfunc(this.srcfile, this.metadata, tspan, ft);

                    % reshape img to X,Y,C,Z,T order
                    img = ImageReader.imreshape(img, info.opts, this.INNER_DIM_ORDER);

                    % write to disk file
                    fwrite(fid, img, "uint16");

                    this.caller.SetProgressBar(k/n_piece);
                end
                fclose(fid);

                if this.state == "off"
                    % remove the exists file
                    delete(this.srcfile);
                    this.metadata = [];
                    this.t = [];
                    this.caller.SetProgressBar(0);
                else
                    this.state = "off";
                    this.data = mpimg(this.folder, [], 1);
                    tmpfile = this.data.FileName;   % 1kb temporary marker
                    this.data.link(dstfile, {"uint16", ...
                        [this.metadata.height, this.metadata.width, ...
                        this.metadata.channels, this.metadata.slices, ...
                        this.metadata.frames], "mov"}, this.INNER_DIM_ORDER);
                    delete(tmpfile);
                end
            else
                % shape as X,Y,C,Z,T
                this.data = zeros([this.metadata.height, this.metadata.width, ...
                        this.metadata.channels, this.metadata.slices, ...
                        this.metadata.frames], "uint16");
                for k = 1:n_piece
                    if this.state == "off"
                        break;
                    end
                    tspan = [(k-1)*block+1, min(k*block, this.metadata.frames)];
                    % load piecewise data
                    img = pfunc(this.srcfile, this.metadata, tspan, ft);

                    % reshape img to X,Y,C,Z,T order
                    img = ImageReader.imreshape(img, info.opts, this.INNER_DIM_ORDER);

                    this.data(:,:,:,:,tspan(1):tspan(2)) = img;
                    this.caller.SetProgressBar(k/n_piece);
                end

                if this.state == "off"
                    this.data = [];
                    this.metadata = [];
                    this.t = [];
                    this.caller.SetProgressBar(0);
                end
            end

            % normalized dimention order
            this.metadata.dimOrder = this.INNER_DIM_ORDER;

            % turn off loader
            this.state = "off";
        end

        function abort(this)
            % This function abort data loading
            this.state = "off";
        end
    end

    methods(Static)
        function tf = fast_loading(ext)
            switch lower(ext)
                case ".nd2"
                    if ispc()
                        try
                            if ~libisloaded('Nd2ReadSdk')
                                [~, ~] = loadlibrary('Nd2ReadSdk', 'Nd2ReadSdk.h');
                            end
                            tf = true;
                        catch
                            tf = false;
                        end
                    else
                        tf = false;
                    end
                case ".ims"
                    tf = true;
                case ".tif"
                    tf = isPyReady();
                otherwise
            end
        end

        function img = imreshape(img, opts, dord)
            frames = numel(img) ...
                /(opts.height*opts.width*opts.channels*opts.slices);

            smap = struct("X", opts.width, ...
                          "Y", opts.height, ...
                          "Z", opts.slices, ...
                          "C", opts.channels, ...
                          "T", frames);

            % reshape
            ss = zeros(1, numel(opts.dimOrder));
            for k = 1:numel(opts.dimOrder), ss(k)=smap.(opts.dimOrder(k)); end
            ss = ss([2,1,3:end]);   % rows as Y, cols as X
            img = reshape(img, ss);

            % permute as "X","Y","C","Z","T"
            [~, new_order] = ismember(opts.dimOrder, dord);
            img = permute(img, new_order);
        end
    end
end

