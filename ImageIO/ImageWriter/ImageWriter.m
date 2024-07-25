classdef ImageWriter < handle
    %IMAGEWRITER This class is an image writer defination, which  support
    %mutiple image format writting
    
    properties(Access=private, Hidden)
        caller
        file
        status  (1,1)   double  {mustBeMember(status, [0, -1])} = 0
        state   (1,1)   string  {mustBeMember(state, ["on","off"])} = "off"
    end

    properties(Access=public, Dependent)
        Status
        State
        Path
    end
    
    methods
        function this = ImageWriter(caller_, file_)
            arguments
                caller_ (1,1)   Register
                file_   (1,1)   string
            end

            this.caller = caller_;
            this.file = file_;
        end

        function r = get.State(this)
            r = this.state;
        end

        function r = get.Status(this)
            r = this.status;
        end

        function r = get.Path(this)
            [r, ~, ~] = fileparts(this.file);
        end
    end

    methods(Access=public)
        function save(this, mov, ch, metadata, block)
            arguments
                this
                mov         (1,1)   regmov
                ch          (1,1)   string  {mustBeMember(ch, ["r","g","b"])}
                metadata    (1,1)   struct
                block       (1,1)   double {mustBePositive, mustBeInteger} = 50
            end
            
            [path, fname, ext] = fileparts(this.file);

            % reselect the save file name
            if ~ismember(ext, [".tif",".hdf5", ".mat"]) || ~isfolder(path)
                % uigetfile seletion
                [dstfile, path] = uiputfile(...
                    {'*.tif','Tag Image File Format(*.tif)'; ...
                    '*.hdf5', 'Hierarchical Data Format(*.hdf5)'; ...
                    '*.mat', 'MATLAB File Format(*.mat)'},...
                    'volumes series selector', fname);
                if isequal(dstfile,0) || isequal(path,0)
                    this.state = "off";
                    return;
                else
                    this.file = string(fullfile(path, dstfile));
                end
            end
            [~, ~, ext] = fileparts(this.file);

            % import the data writter setting
            pathToWriter = fileparts(mfilename('fullpath'));
            rf = importdata([pathToWriter, filesep, 'configuration.ini']);
            rf = string(rf).split(":");
            if ~ismember(upper(ext), rf(:,1))
                throw(MException("ImageLoader:invalidImageWrittingFunc", ...
                        "Can not parse function without registration."));
            end

            % pfunc: [meta, data] = pfunc(file, tspan)
            pfunc = str2func(rf(upper(ext)==rf(:,1), 2));

            % if turbo loading enabled
            ft = ImageWriter.fast_writting(ext);

            switch lower(ext)
                case ".tif"
                    % reset the block size
                    block = mov.MetaData.frames;
                    % indeterminate progress bar
                    this.caller.SetProgressBar(0, true);
                otherwise
                    block = min(ceil(mov.MetaData.frames/20), block);
                    this.caller.SetProgressBar(0);
            end

            % get data piecewiselly
            n_piece = ceil(mov.MetaData.frames/block);

            % turn on saver
            this.state = "on";

            for k = 1:n_piece
                if this.state == "off"
                    break;
                end
                tspan = [(k-1)*block+1, min(k*block, mov.MetaData.frames)];

                % save piecewise data
                this.status = pfunc(this.file, mov, ch, tspan, metadata, ft);
                if this.status == -1
                    this.state = "off";
                    break;
                end

                this.caller.SetProgressBar(k/n_piece);
            end

            % turn off writer
            this.state = "off";
        end
    end

    methods(Static)
        function tf = fast_writting(ext)
            switch lower(ext)
                case ".tif"
                    tf = isPyReady();
                otherwise
                    tf = false;
            end
        end
    end
end

