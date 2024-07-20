function info = loadmeta_bf(file)
%LOADMETA_BF This function load the metadata by using package bfmatlab
% Input:
%   - file: 1-by-1 string or char array, image file path
% Output:
%   - info: 1-by-1 struct, with 1-by-12 metadata table and n-by-1
%           real time array, field {opts, rt}

arguments
    file    (1,:)   {mustBeFile}
end

[~, ~, ext] = fileparts(file);

status = bfCheckJavaPath(1);
assert(status, ['Missing Bio-Formats library. Either add bioformats_package.jar '...
    'to the static Java path or add it to the Matlab path.']);
bfInitLogging();
r = bfGetReader(char(file), 0);

info_bf = r.getMetadataStore();

% extract the information
width = info_bf.getPixelsSizeX(0).getValue();
height = info_bf.getPixelsSizeY(0).getValue();
channels = info_bf.getPixelsSizeC(0).getValue();
slices = info_bf.getPixelsSizeZ(0).getValue();
frames = info_bf.getPixelsSizeT(0).getValue();
images = channels*slices*frames;
dimOrder = string(split(info_bf.getPixelsDimensionOrder(0).getValue(),""));
dimOrder = dimOrder(strlength(dimOrder)>0)';
cOrder = get_channel_order(info_bf, ext);
dataType = string(info_bf.getPixelsType(0).getValue());
if dataType == "float"
    dataType = "single";
end

switch lower(ext)
    case ".tif"
        % Ask the user to input necessary information
        prompt = {'Enter the z scan thickness(\mum):',...
            'Enter the binning size:',...
            'Enter total optical magnification:',...
            'Enter the camera pixel size(\mum):'};
        dlgtitle = 'Microscope parameters setting';
        dims = [1 40];
        definput = {'1.5';'2';'60';'6.5'};
        dlgopts.Interpreter = 'tex';
        dlgout = inputdlg(prompt,dlgtitle,dims,definput,dlgopts);
        if isempty(dlgout), dlgout = definput; end

        % extract and calculate the X,Y resolution
        zRes = str2double(dlgout{1});
        binsz = str2double(dlgout{2});
        tcom = str2double(dlgout{3});
        cps = str2double(dlgout{4});

        % since the camera physical pixel is square,
        % and the binning size is also square
        xRes = cps/tcom*binsz;
        yRes = xRes;
    case ".ims"
        xRes = double(info_bf.getPixelsPhysicalSizeX(0).value);
        yRes = double(info_bf.getPixelsPhysicalSizeY(0).value);
        if slices ~= 1
            zRes = double(info_bf.getPixelsPhysicalSizeZ(0).value) ...
                *(slices/(slices-1));
        else
            zRes = 0;   % only one plane, no resolution at z direction
        end
    case ".nd2"
        xRes = double(info_bf.getPixelsPhysicalSizeX(0).value);
        yRes = double(info_bf.getPixelsPhysicalSizeY(0).value);
        if slices ~= 1
            zRes = double(info_bf.getPixelsPhysicalSizeZ(0).value);
        else
            zRes = 0;   % only one plane, no resolution at z direction
        end
    otherwise

end

opts = table(width,...          % #pixel
            height,...          % #pixel
            channels,...        % #channels
            slices,...          % #slices
            frames,...          % #volumes
            images,...          % #images
            xRes,...            % μm/pixel
            yRes,...            % μm/pixel
            zRes,...            % μm/layer z
            dataType,...        % uint8, uint16, single
            dimOrder,...        % dimension order array
            cOrder);            % color channel order

rt = get_real_time(r, ext);

r.close();

info = struct("opts", opts, "rt", rt);
end

function r = get_channel_order(info_bf, ext)
switch lower(ext)
    case ".tif"
        cc = info_bf.getChannelCount(0);
        C = zeros(cc,4,'uint8');
        for k = 0:cc-1
            cinfo = info_bf.getChannelColor(0,k);
            % some format there is no color channel information such as
            % tiff file exported from Fiji
            if isempty(cinfo)
                disp("The channel order not found. Manual inputs are required.");
                cn = info_bf.getPixelsSizeC(0).getValue(); % color channels number
                cdlgtitle = 'Camera color channel setting';
                cdims = [1 40];
                switch cn
                    case 1
                        cprompt = {'Enter the channel color(r/g/b):'};
                        cdefinput = {'g'};
                    case 2
                        cprompt = {'Enter the channel-1 color(r/g/b):',...
                            'Enter the channel-2 color(r/g/b):'};
                        cdefinput = {'r','g'};
                    otherwise
                        error("Unsupported Color Channel.");
                end
                cdlgout = inputdlg(cprompt,cdlgtitle,cdims,cdefinput);
                if isempty(cdlgout) || numel(cdlgout) ~= cn
                    cdlgout = cdefinput;
                end
                r = reshape(string(cdlgout),1,[]);

                return;
            end
            C(k+1,:) = [cinfo.getRed(),cinfo.getGreen(),cinfo.getBlue(),...
                cinfo.getAlpha()];
        end
        C_na = C(:,1:3);  % cut the color map without alpha channel
        % select the maximum value to represent the channel
        [~,order] = max(C_na,[],2);
        r = ["r","g","b"];
        r = r(order);
    otherwise
        cc = info_bf.getChannelCount(0);
        C = zeros(cc,4,'uint8');
        for k = 0:cc-1
            cinfo = info_bf.getChannelColor(0,k);
            % some format there is no color channel information such as
            % tiff file exported from Fiji
            if isempty(cinfo)
                r = "";
                return;
            end
            C(k+1,:) = [cinfo.getRed(),cinfo.getGreen(),cinfo.getBlue(),...
                cinfo.getAlpha()];
        end
        C_na = C(:,1:3);  % cut the color map without alpha channel
        % select the maximum value to represent the channel
        [~,order] = max(C_na,[],2);
        r = ["r","g","b"];
        r = r(order);
end

end

function t = get_real_time(r, ext)
switch lower(ext)
    case ".tif"
        % get packaged meta data
        gmd = r.getGlobalMetadata();
        smd = r.getMetadataStore();
        frames = smd.getPixelsSizeT(0).getValue();

        % get delta time relative to world time
        if ~isempty(gmd.get("AcquisitionStart")) ...
                && ~isempty(gmd.get("RecordingDate"))
            dt = seconds(diff([datetime(gmd.get("AcquisitionStart")),...
                datetime(gmd.get("RecordingDate"))]));
        else
            warning("The acquire time not found. Frame Indices replaced.");
            t = [];
            return;
        end

        if isa(gmd,'java.util.Hashtable')
            keys = gmd.keys();
            n = 1;
            t = [];
            while keys.hasMoreElements()
                tp_name = ['TimePoint',num2str(n)];
                if gmd.containsKey(tp_name)
                    t = [t,datetime(gmd.get(tp_name))]; %#ok<AGROW>
                    n = n + 1;
                end
                keys.nextElement();    % inner iterator increase
            end
        else
            error('invalid input format');
        end

        % datetime array t, change it to seconds
        t = [0,cumsum(seconds(diff(t)))]';
        t(frames+1:end) = [];
        if frames > 1
            t(1) = (t(end)-t(2))/(frames - 2 + eps);   % average estimation
            t = t + dt;
        end

        % control the precision: ms
        t = round(t, 3);
    case ".nd2"
        % get packaged meta data
        smd = r.getMetadataStore();
        channels = smd.getPixelsSizeC(0).getValue();
        slices = smd.getPixelsSizeZ(0).getValue();
        frames = smd.getPixelsSizeT(0).getValue();

        t = nan(frames, 1);
        for k = 1:frames-1
            t(k+1) = smd.getPlaneDeltaT(0,slices*channels*k-1).value();
        end
        t = fillmissing(t, "linear");
        t = round(t, 3);
    case ".ims"
        % get packaged meta data
        gmd = r.getGlobalMetadata();
        smd = r.getMetadataStore();
        frames = smd.getPixelsSizeT(0).getValue();

        t0 = datetime(gmd.get("TimePoint1"));
        t = nan(frames, 1);
        for k = 1:frames
            t(k) = seconds(datetime(gmd.get(sprintf("TimePoint%d", k))) ...
                - t0);
        end
        
        t = round(t, 3);
    otherwise
end
end