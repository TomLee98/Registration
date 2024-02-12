function [info, data] = loadims(file, tspan)
%LOADIMS This function load ims files with tspan
% Input:
%   - file: 1-by-1 string, the data full file path
%   - tspan: 1-by-2 positive integer array, [start, end] of frame indices
% Output:
%   - info: struct with {opts, rt}, where opts is 1-by-12 metadata table,
%     rt is n-by-1 nonnegtive double time array
%
% if date is not needed but info, use
%       info = loadims(file)
% or you want the data, use
%       [info, data] = loadims(file, tspan)
% where you can pass the time span as tspan for continous time cut

arguments
    file    (1,1)   string
    tspan   (1,:)   double = []
end

% generate HDF5 reader
GID = H5G.open(H5F.open(file), '/DataSet');
dataobj = DatasetReader(GID);

nargoutchk(1,2);

% extract the information
width = dataobj.SizeX;
height = dataobj.SizeY;
channels = dataobj.SizeC;
slices = dataobj.SizeZ;
frames  = dataobj.SizeT;
images = channels*slices*frames;
xRes = (dataobj.ExtendMaxX-dataobj.ExtendMinX)/dataobj.SizeX;   % um/pix
yRes = (dataobj.ExtendMaxY-dataobj.ExtendMinY)/dataobj.SizeY;   % um/pix
zRes = (dataobj.ExtendMaxZ-dataobj.ExtendMinZ)/(dataobj.SizeZ-1);   % um/pix
dataType = string(dataobj.DataType);
if dataType == "float"
    dataType = "single";
end
dimOrder = ["X","Y","Z","C","T"];
cOrder = getChannelOrder(dataobj.ChannelInfo);

% generate opts
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
            dimOrder,...        % dimention order array
            cOrder);            % color channel order

% extract time stamps
rt = getRelativeTime(dataobj.Timestamps);

info.opts = opts;
info.rt = rt;

if nargout == 2
    % load the image data
    data = dataobj.GetData(tspan);
end

    function cOrder = getChannelOrder(cInfo)
        order = ["r","g","b"];
        cOrder = strings(1, numel(cInfo));
        for k = 1:numel(cInfo)
            [~, p] = max(cInfo(k).Color);
            cOrder(k) = order(p);
        end
    end

    function rt = getRelativeTime(timeStamp)
        rt = datetime(string(timeStamp),"Format","uuuu-MM-dd HH:mm:ss.SSS");
        rt = seconds(rt - rt(1));
    end
end

